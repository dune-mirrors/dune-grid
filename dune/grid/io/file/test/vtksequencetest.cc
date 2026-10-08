// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <cmath>
#include <condition_variable>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <regex>
#include <string>
#include <thread>
#include <vector>

#if __has_include(<execinfo.h>) && __has_include(<pthread.h>) && __has_include(<unistd.h>)
#define VTKSEQUENCETEST_HAVE_WATCHDOG 1
#include <execinfo.h>
#include <pthread.h>
#include <unistd.h>
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/path.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

// The test occasionally timed out in CI without any hint where it got stuck.
// To diagnose this, we keep a description of the current step, and print it
// together with a backtrace of the main thread when the test receives SIGUSR1,
// either sent by ctest on timeout (TIMEOUT_SIGNAL_NAME, CMake >= 3.27) or by
// the watchdog thread below.
namespace {

char currentStep[256] = "initialization";

void setStep(int dim, const std::string& name, int step)
{
  std::snprintf(currentStep, sizeof(currentStep), "dim=%d, sequence=%s, step=%d", dim, name.c_str(), step);
}

#if VTKSEQUENCETEST_HAVE_WATCHDOG
extern "C" void dumpStateAndExit(int)
{
  static const char msg[] = "\nvtksequencetest got stuck at: ";
  [[maybe_unused]] auto r1 = ::write(STDERR_FILENO, msg, sizeof(msg)-1);
  [[maybe_unused]] auto r2 = ::write(STDERR_FILENO, currentStep, strnlen(currentStep, sizeof(currentStep)));
  [[maybe_unused]] auto r3 = ::write(STDERR_FILENO, "\n", 1);
  void* frames[64];
  int n = backtrace(frames, 64);
  backtrace_symbols_fd(frames, n, STDERR_FILENO);
  _exit(1);
}

// Sends SIGUSR1 to the main thread if the test does not finish in time
class Watchdog
{
public:
  explicit Watchdog(std::chrono::seconds timeout)
    : mainThread_(pthread_self())
  {
    // make sure that libgcc is loaded before backtrace() is called inside the signal handler
    void* frame;
    backtrace(&frame, 1);
    std::signal(SIGUSR1, dumpStateAndExit);

    thread_ = std::thread([this,timeout]{
      std::unique_lock lock(mutex_);
      if (!cv_.wait_for(lock, timeout, [this]{ return done_; }))
        pthread_kill(mainThread_, SIGUSR1);
    });
  }

  ~Watchdog()
  {
    {
      std::lock_guard lock(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    thread_.join();
  }

private:
  pthread_t mainThread_;
  std::thread thread_;
  std::mutex mutex_;
  std::condition_variable cv_;
  bool done_ = false;
};
#endif

} // end anonymous namespace


// A time-dependent vertex function, to have data that changes between the time steps
template<class GridView>
class TimeDependentFunction
  : public Dune::VTKWriter<GridView>::VTKFunction
{
  using Entity = typename GridView::template Codim<0>::Entity;
  using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

public:
  void setTime (double time) { time_ = time; }

  int ncomps () const override { return 1; }

  double evaluate ([[maybe_unused]] int comp, [[maybe_unused]] const Entity& e,
                   [[maybe_unused]] const LocalCoordinate& xi) const override
  {
    return std::sin(time_);
  }

  std::string name () const override { return "timeFunction"; }

private:
  double time_ = 0.0;
};


struct DataSet
{
  double timestep;
  std::string file;
};

// Extract the `timestep` and `file` attributes of all DataSet entries of a pvd file
std::vector<DataSet> readPvdFile (const std::string& filename)
{
  std::ifstream pvdFile(filename);
  if (!pvdFile)
    DUNE_THROW(Dune::IOError, "File " << filename << " could not be opened!");

  static const std::regex dataSetRegex(R"re(<DataSet\s+timestep="([^"]*)".*\sfile="([^"]*)")re");

  std::vector<DataSet> dataSets;
  std::string line;
  std::smatch match;
  while (std::getline(pvdFile, line))
    if (std::regex_search(line, match, dataSetRegex))
      dataSets.push_back({std::stod(match[1]), match[2]});
  return dataSets;
}

// Check that the pvd file lists the given time steps and that all referenced files exist
template<class GridView>
Dune::TestSuite checkPvdFile (const std::string& name, const std::string& path,
                              const std::vector<double>& timesteps)
{
  Dune::TestSuite t("checkPvdFile(" + name + ")");

  const std::string extension = GridView::dimension == 1 ? ".vtp" : ".vtu";
  const auto dataSets = readPvdFile(name + ".pvd");
  t.require(dataSets.size() == timesteps.size())
    << name << ".pvd contains " << dataSets.size() << " data sets, expected " << timesteps.size();

  for (std::size_t i = 0; i < dataSets.size(); ++i) {
    t.check(std::abs(dataSets[i].timestep - timesteps[i]) < 1e-8)
      << name << ".pvd: wrong timestep " << dataSets[i].timestep << " of data set " << i
      << ", expected " << timesteps[i];

    char seqName[32];
    std::snprintf(seqName, sizeof(seqName), "-%05zu", i);
    const auto expectedFile = Dune::concatPaths(path, name + seqName + extension);
    t.check(dataSets[i].file == expectedFile)
      << name << ".pvd: data set " << i << " references " << dataSets[i].file
      << ", expected " << expectedFile;
    t.check(std::filesystem::exists(dataSets[i].file))
      << name << ".pvd: referenced file " << dataSets[i].file << " does not exist";
  }

  return t;
}

// Write a sequence of `numSteps` time steps. If a `restartWriter` is given,
// continue its time sequence, as is done when restarting a simulation.
template<class GridView>
std::vector<double> writeSequence (const GridView& gridView, const std::string& name,
                                   const std::string& path, int numSteps,
                                   Dune::VTK::OutputType type = Dune::VTK::ascii,
                                   const std::vector<double>& previousTimesteps = {})
{
  constexpr int dim = GridView::dimension;
  std::vector<int> vertexData(gridView.size(dim), dim);
  std::vector<int> cellData(gridView.size(0), 0);
  auto timeFunction = std::make_shared<TimeDependentFunction<GridView>>();

  auto vtkWriter = std::make_shared<Dune::VTKWriter<GridView>>(gridView);
  Dune::VTKSequenceWriter<GridView> vtk(vtkWriter, name, path, "");
  vtk.setTimeSteps(previousTimesteps);
  vtk.addVertexData(vertexData, "vertexData");
  vtk.addCellData(cellData, "cellData");
  vtk.addVertexData(timeFunction);

  double time = previousTimesteps.empty() ? 0.0 : previousTimesteps.back();
  for (int i = 0; i < numSteps; ++i) {
    time += 0.25;
    setStep(dim, name, i);
    timeFunction->setTime(time);
    vtk.write(time, type);
  }

  return vtk.getTimeSteps();
}

template<int dim>
Dune::TestSuite vtkSequenceCheck ()
{
  std::cout << "vtkSequenceCheck dim=" << dim << std::endl;
  Dune::TestSuite t("vtkSequenceCheck<" + std::to_string(dim) + ">");

  Dune::FieldVector<double,dim> h(1.0);
  std::array<int,dim> n;
  n.fill(2);
  Dune::YaspGrid<dim> grid(h, n);
  const auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);

  const std::string prefix = "vtksequencetest-" + std::to_string(dim) + "D";

  { // simple sequence in the current directory
    const auto name = prefix + "-ascii";
    const auto timesteps = writeSequence(gridView, name, ".", 3);
    t.check(timesteps == std::vector<double>{0.25, 0.5, 0.75}) << "wrong time steps returned by getTimeSteps()";
    t.subTest(checkPvdFile<GridView>(name, ".", timesteps));
  }

  { // other output type and pieces written to a subdirectory
    const auto name = prefix + "-subdir";
    const std::string path = "vtksequencetest-output";
    std::filesystem::create_directories(path);
    const auto timesteps = writeSequence(gridView, name, path, 2, Dune::VTK::appendedraw);
    t.subTest(checkPvdFile<GridView>(name, path, timesteps));
  }

  { // restart: a second writer continues the sequence of the first one
    const auto name = prefix + "-restart";
    const auto timesteps1 = writeSequence(gridView, name, ".", 2);
    const auto timesteps2 = writeSequence(gridView, name, ".", 2, Dune::VTK::ascii, timesteps1);
    t.check(timesteps2 == std::vector<double>{0.25, 0.5, 0.75, 1.0}) << "wrong time steps after restart";
    t.subTest(checkPvdFile<GridView>(name, ".", timesteps2));
  }

  return t;
}

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

#if VTKSEQUENCETEST_HAVE_WATCHDOG
  Watchdog watchdog(std::chrono::seconds(120));
#endif

  Dune::TestSuite t;
  t.subTest(vtkSequenceCheck<1>());
  t.subTest(vtkSequenceCheck<2>());
  t.subTest(vtkSequenceCheck<3>());
  setStep(0, "finalization", 0);
  return t.exit();
}
