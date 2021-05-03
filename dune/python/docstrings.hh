/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by pybind11_mkdoc.
 */

#define __EXPAND(x)                                      x
#define __COUNT(_1, _2, _3, _4, _5, _6, _7, COUNT, ...)  COUNT
#define __VA_SIZE(...)                                   __EXPAND(__COUNT(__VA_ARGS__, 7, 6, 5, 4, 3, 2, 1))
#define __CAT1(a, b)                                     a ## b
#define __CAT2(a, b)                                     __CAT1(a, b)
#define __DOC1(n1)                                       __doc_##n1
#define __DOC2(n1, n2)                                   __doc_##n1##_##n2
#define __DOC3(n1, n2, n3)                               __doc_##n1##_##n2##_##n3
#define __DOC4(n1, n2, n3, n4)                           __doc_##n1##_##n2##_##n3##_##n4
#define __DOC5(n1, n2, n3, n4, n5)                       __doc_##n1##_##n2##_##n3##_##n4##_##n5
#define __DOC6(n1, n2, n3, n4, n5, n6)                   __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6
#define __DOC7(n1, n2, n3, n4, n5, n6, n7)               __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6##_##n7
#define DOC(...)                                         __EXPAND(__EXPAND(__CAT2(__DOC, __VA_SIZE(__VA_ARGS__)))(__VA_ARGS__))

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif


static const char *__doc_Dune_MultipleCodimMultipleGeomTypeMapper = R"doc()doc";

static const char *__doc_Dune_Python_BoundarySegment = R"doc()doc";

static const char *__doc_Dune_Python_BoundarySegment_BoundarySegment = R"doc()doc";

static const char *__doc_Dune_Python_BoundarySegment_operator_call = R"doc()doc";

static const char *__doc_Dune_Python_BoundarySegment_parametrization = R"doc()doc";

static const char *__doc_Dune_Python_Capabilities_HasGridFactory = R"doc()doc";

static const char *__doc_Dune_Python_Capabilities_HasStructuredGridFactory = R"doc()doc";

static const char *__doc_Dune_Python_Capabilities_canIterate = R"doc()doc";

static const char *__doc_Dune_Python_EvaluateType = R"doc()doc";

static const char *__doc_Dune_Python_EvaluateType_name = R"doc()doc";

static const char *__doc_Dune_Python_FunctionRange = R"doc()doc";

static const char *__doc_Dune_Python_FunctionRange_value = R"doc()doc";

static const char *__doc_Dune_Python_GridFunctionTraits = R"doc()doc";

static const char *__doc_Dune_Python_GridFunctionTraits_2 = R"doc()doc";

static const char *__doc_Dune_Python_GridFunctionTraits_3 = R"doc()doc";

static const char *__doc_Dune_Python_GridModificationListener = R"doc()doc";

static const char *__doc_Dune_Python_GridModificationListener_postModification = R"doc()doc";

static const char *__doc_Dune_Python_GridModificationListener_preModification = R"doc()doc";

static const char *__doc_Dune_Python_GridObjectTraits = R"doc()doc";

static const char *__doc_Dune_Python_GridViewPartition = R"doc()doc";

static const char *__doc_Dune_Python_GridViewPartition_GridViewPartition = R"doc()doc";

static const char *__doc_Dune_Python_GridViewPartition_gridView = R"doc()doc";

static const char *__doc_Dune_Python_GridViewPartition_obj = R"doc()doc";

static const char *__doc_Dune_Python_LocalEvaluatorAdapter = R"doc()doc";

static const char *__doc_Dune_Python_LocalEvaluatorAdapter_LocalEvaluatorAdapter = R"doc()doc";

static const char *__doc_Dune_Python_LocalEvaluatorAdapter_evaluator = R"doc()doc";

static const char *__doc_Dune_Python_LocalEvaluatorAdapter_operator_call = R"doc()doc";

static const char *__doc_Dune_Python_LocalEvaluatorAdapter_operator_call_2 = R"doc()doc";

static const char *__doc_Dune_Python_Marker = R"doc()doc";

static const char *__doc_Dune_Python_Marker_Coarsen = R"doc()doc";

static const char *__doc_Dune_Python_Marker_Keep = R"doc()doc";

static const char *__doc_Dune_Python_Marker_Refine = R"doc()doc";

static const char *__doc_Dune_Python_NumPyCommDataHandle = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_ProxyDataHandle = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_contains = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_contains_2 = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_fixedSize = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_fixedSze = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_gather = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_gather_2 = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_method = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_scatter = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_scatter_2 = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_size = R"doc()doc";

static const char *__doc_Dune_Python_ProxyDataHandle_size_2 = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_PyBoundaryIntersectionIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_elementIt = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_gridView = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_intersectionEnd = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_intersectionIt = R"doc()doc";

static const char *__doc_Dune_Python_PyBoundaryIntersectionIterator_next = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_PyGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_PyGridFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_bind = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_evaluate = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_lf = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_operator_call = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_pyObj = R"doc()doc";

static const char *__doc_Dune_Python_PyGridFunction_unbind = R"doc()doc";

static const char *__doc_Dune_Python_PyHierarchicIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyHierarchicIterator_PyHierarchicIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyHierarchicIterator_end = R"doc()doc";

static const char *__doc_Dune_Python_PyHierarchicIterator_it = R"doc()doc";

static const char *__doc_Dune_Python_PyHierarchicIterator_next = R"doc()doc";

static const char *__doc_Dune_Python_PyIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyIterator_PyIterator = R"doc()doc";

static const char *__doc_Dune_Python_PyIterator_end = R"doc()doc";

static const char *__doc_Dune_Python_PyIterator_it = R"doc()doc";

static const char *__doc_Dune_Python_PyIterator_next = R"doc()doc";

static const char *__doc_Dune_Python_Reader = R"doc()doc";

static const char *__doc_Dune_Python_Reader_dgf = R"doc()doc";

static const char *__doc_Dune_Python_Reader_dgfString = R"doc()doc";

static const char *__doc_Dune_Python_Reader_gmsh = R"doc()doc";

static const char *__doc_Dune_Python_Reader_structured = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGlobalGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGlobalGridFunction_SimpleGlobalGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGlobalGridFunction_operator_call = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction_SimpleGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction_gridView = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction_gridView_2 = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction_localEvaluator = R"doc()doc";

static const char *__doc_Dune_Python_SimpleGridFunction_localEvaluator_2 = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_SimpleLocalFunction = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_SimpleLocalFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_bind = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_element = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_element_2 = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_localEvaluator = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_operator_call = R"doc()doc";

static const char *__doc_Dune_Python_SimpleLocalFunction_unbind = R"doc()doc";

static const char *__doc_Dune_Python_VTKDataType = R"doc()doc";

static const char *__doc_Dune_Python_VTKDataType_CellData = R"doc()doc";

static const char *__doc_Dune_Python_VTKDataType_CellVector = R"doc()doc";

static const char *__doc_Dune_Python_VTKDataType_PointData = R"doc()doc";

static const char *__doc_Dune_Python_VTKDataType_PointVector = R"doc()doc";

static const char *__doc_Dune_Python_addToVTKWriter = R"doc()doc";

static const char *__doc_Dune_Python_coordinates = R"doc()doc";

static const char *__doc_Dune_Python_coordinates_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_CommOp = R"doc()doc";

static const char *__doc_Dune_Python_detail_CommOp_add = R"doc()doc";

static const char *__doc_Dune_Python_detail_CommOp_set = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertBoundaries = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertBoundaries_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_4 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_5 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_6 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertElements_7 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertVertices = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertVertices_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertVertices_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridFactory_insertVertices_4 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridId = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridId_id = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridObjectTraits_element = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridObjectTraits_element_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridObjectTraits_localCoordinate = R"doc()doc";

static const char *__doc_Dune_Python_detail_GridObjectTraits_localCoordinate_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_LocalViewRegistry = R"doc()doc";

static const char *__doc_Dune_Python_detail_LocalViewRegistry_bind = R"doc()doc";

static const char *__doc_Dune_Python_detail_LocalViewRegistry_binds = R"doc()doc";

static const char *__doc_Dune_Python_detail_LocalViewRegistry_find = R"doc()doc";

static const char *__doc_Dune_Python_detail_LocalViewRegistry_unbind = R"doc()doc";

static const char *__doc_Dune_Python_detail_PyGridFunctionEvaluator = R"doc()doc";

static const char *__doc_Dune_Python_detail_addGridModificationListener = R"doc()doc";

static const char *__doc_Dune_Python_detail_callLocalFunction = R"doc()doc";

static const char *__doc_Dune_Python_detail_callLocalFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_callLocalFunction_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_gridModificationListeners = R"doc()doc";

static const char *__doc_Dune_Python_detail_gridModificationListeners_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_gridView = R"doc()doc";

static const char *__doc_Dune_Python_detail_gridView_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_gridView_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_idSetSubId = R"doc()doc";

static const char *__doc_Dune_Python_detail_indexSetSubIndex = R"doc()doc";

static const char *__doc_Dune_Python_detail_localViewRegistry = R"doc()doc";

static const char *__doc_Dune_Python_detail_makeSubEntities = R"doc()doc";

static const char *__doc_Dune_Python_detail_makeSubEntity = R"doc()doc";

static const char *__doc_Dune_Python_detail_mapperCommunicate = R"doc()doc";

static const char *__doc_Dune_Python_detail_mapperSubIndex = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_eq = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_ge = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_gt = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_le = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_le_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_lt = R"doc()doc";

static const char *__doc_Dune_Python_detail_operator_ne = R"doc()doc";

static const char *__doc_Dune_Python_detail_pushForwardGradients = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerExtendedEntityInterface = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerExtendedEntityInterface_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerExtendedEntityInterface_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridEntity = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridGeometry = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridIntersection = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridViewConstructorFromGrid = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridViewConstructorFromGrid_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerGridViewConstructorFromGrid_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerMapperSubIndex = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerMapperSubIndex_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerMapperSubIndex_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerPyGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubId = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubId_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubId_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubIndex = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubIndex_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_registerSubIndex_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_to_string = R"doc()doc";

static const char *__doc_Dune_Python_detail_to_string_2 = R"doc()doc";

static const char *__doc_Dune_Python_detail_to_string_3 = R"doc()doc";

static const char *__doc_Dune_Python_detail_to_string_4 = R"doc()doc";

static const char *__doc_Dune_Python_fillGridFactory = R"doc()doc";

static const char *__doc_Dune_Python_flatCopy = R"doc()doc";

static const char *__doc_Dune_Python_flatCopy_2 = R"doc()doc";

static const char *__doc_Dune_Python_gridView = R"doc()doc";

static const char *__doc_Dune_Python_makeMultipleCodimMultipleGeomTypeMapper = R"doc()doc";

static const char *__doc_Dune_Python_makeNumPyArray = R"doc()doc";

static const char *__doc_Dune_Python_makePyGridViewIterator = R"doc()doc";

static const char *__doc_Dune_Python_makePyGridViewIterator_2 = R"doc()doc";

static const char *__doc_Dune_Python_makePyGridViewPartitionIterator = R"doc()doc";

static const char *__doc_Dune_Python_makePyGridViewPartitionIterator_2 = R"doc()doc";

static const char *__doc_Dune_Python_numPyCommDataHandle = R"doc()doc";

static const char *__doc_Dune_Python_numPyCommDataHandle_2 = R"doc()doc";

static const char *__doc_Dune_Python_pyGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_pyGridFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_readDGF = R"doc()doc";

static const char *__doc_Dune_Python_readDGF_2 = R"doc()doc";

static const char *__doc_Dune_Python_readGmsh = R"doc()doc";

static const char *__doc_Dune_Python_readGmsh_2 = R"doc()doc";

static const char *__doc_Dune_Python_reader = R"doc()doc";

static const char *__doc_Dune_Python_reader_2 = R"doc()doc";

static const char *__doc_Dune_Python_reader_3 = R"doc()doc";

static const char *__doc_Dune_Python_reader_4 = R"doc()doc";

static const char *__doc_Dune_Python_reader_5 = R"doc()doc";

static const char *__doc_Dune_Python_reader_6 = R"doc()doc";

static const char *__doc_Dune_Python_reader_7 = R"doc()doc";

static const char *__doc_Dune_Python_registerDataHandle = R"doc()doc";

static const char *__doc_Dune_Python_registerGridEntities = R"doc()doc";

static const char *__doc_Dune_Python_registerGridEntities_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerGridEntity = R"doc()doc";

static const char *__doc_Dune_Python_registerGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_registerGridFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerGridFunction_3 = R"doc()doc";

static const char *__doc_Dune_Python_registerGridGeometry = R"doc()doc";

static const char *__doc_Dune_Python_registerGridIdSet = R"doc()doc";

static const char *__doc_Dune_Python_registerGridIntersection = R"doc()doc";

static const char *__doc_Dune_Python_registerGridView = R"doc()doc";

static const char *__doc_Dune_Python_registerGridViewIndexSet = R"doc()doc";

static const char *__doc_Dune_Python_registerGridViewIndexSet_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerGridViewPartition = R"doc()doc";

static const char *__doc_Dune_Python_registerHierarchicalGrid = R"doc()doc";

static const char *__doc_Dune_Python_registerHierarchicalGridIdSets = R"doc()doc";

static const char *__doc_Dune_Python_registerHierarchicalGridPicklingSupport = R"doc()doc";

static const char *__doc_Dune_Python_registerHierarchicalGridPicklingSupport_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerHierarchicalGridPicklingSupport_3 = R"doc()doc";

static const char *__doc_Dune_Python_registerLocalView = R"doc()doc";

static const char *__doc_Dune_Python_registerMapper = R"doc()doc";

static const char *__doc_Dune_Python_registerMapperCommunicate = R"doc()doc";

static const char *__doc_Dune_Python_registerMultipleCodimMultipleGeomTypeMapper = R"doc()doc";

static const char *__doc_Dune_Python_registerMultipleCodimMultipleGeomTypeMapper_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerPersistentContainer = R"doc()doc";

static const char *__doc_Dune_Python_registerPyBoundaryIntersectionIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewIterator_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewIterator_3 = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewPartitionIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewPartitionIterator_2 = R"doc()doc";

static const char *__doc_Dune_Python_registerPyGridViewPartitionIterator_3 = R"doc()doc";

static const char *__doc_Dune_Python_registerPyHierarchicIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerPyIntersectionIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerPyIterator = R"doc()doc";

static const char *__doc_Dune_Python_registerVTKWriter = R"doc()doc";

static const char *__doc_Dune_Python_simpleGridFunction = R"doc()doc";

static const char *__doc_Dune_Python_simpleGridFunction_2 = R"doc()doc";

static const char *__doc_Dune_Python_simpleGridFunction_3 = R"doc()doc";

static const char *__doc_Dune_Python_stdFunction = R"doc()doc";

static const char *__doc_Dune_SingletonStorage = R"doc()doc";

static const char *__doc_Dune_SingletonStorage_Item = R"doc()doc";

static const char *__doc_Dune_SingletonStorage_ItemWrapper = R"doc()doc";

static const char *__doc_Dune_SingletonStorage_ItemWrapper_ItemWrapper = R"doc()doc";

static const char *__doc_Dune_SingletonStorage_ItemWrapper_obj = R"doc()doc";

static const char *__doc_Dune_SingletonStorage_instance = R"doc(return singleton instance of given Object type.)doc";

static const char *__doc_Dune_SingletonStorage_storage = R"doc()doc";

static const char *__doc_Dune_YaspGrid = R"doc()doc";

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif
