add_executable(
  main
  BaseChannel.cpp
  BaseDecoder.cpp
  BaseSimulator.cpp
  BpskAwgnChannel.cpp
  BpskBscChannel.cpp
  CommonTransformations.cpp
  ConfigReading.cpp
  CRC.cpp
  DecoderType.cpp
  Encoder.cpp
  GaussianApproximation.cpp
  main.cpp
  MonteCarloSimulator.cpp
  PolarCode.cpp
  ScCrcAidedDecoder.cpp
  ScCreeperDecoder.cpp
  ScDecoder.cpp
  ScDecoderTreeMaker.cpp
  ScFanoDecoder.cpp
  ScFlipDecoder.cpp
  ScFlipFanoDecoder.cpp
  ScFlipProgDecoder.cpp
  ScListDecoder.cpp
  ScListFanoDecoder.cpp
  ScListFlipOracleStatDecoder.cpp
  ScListFlipStatDecoder.cpp
  ScOptimized.cpp
  ScRecursiveDecoder.cpp
  ScStackDecoder.cpp
  ScStackOperationCounter.cpp
  ScStackOptimizedDecoder.cpp
  Simulate.cpp
  SimulationType.cpp
  Test_AwgnBer.cpp
  Test_CRC.cpp
  Test_GA.cpp
  Test_ScTreeMaker.cpp)
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(main PUBLIC default) #mexbind0x::mexbind0x)

if(Matlab_FOUND)
  matlab_add_mex(
    NAME polar_decoder R2017b
    SRC ScCreeperDecoder.cpp
        mex.cpp
        ScOptimized.cpp
        ScStackOptimizedDecoder.cpp
        ScFlipFanoDecoder.cpp
        ScFanoDecoder.cpp
        ScDecoder.cpp
        ScCrcAidedDecoder.cpp
        ScStackDecoder.cpp
        BaseDecoder.cpp
        PolarCode.cpp
        GaussianApproximation.cpp
        CommonTransformations.cpp
        CRC.cpp
    LINK_TO mexbind0x::mexbind0x default)
  target_include_directories(polar_decoder
                             PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
endif()

add_executable(
  test_creeper
  test_creeper.cpp
    ScCreeperDecoder.cpp
        ScOptimized.cpp
        ScStackOptimizedDecoder.cpp
        ScFlipFanoDecoder.cpp
        ScFanoDecoder.cpp
        ScDecoder.cpp
        ScCrcAidedDecoder.cpp
        ScStackDecoder.cpp
        BaseDecoder.cpp
        PolarCode.cpp
        GaussianApproximation.cpp
        CommonTransformations.cpp
        CRC.cpp
        )

target_link_libraries(test_creeper PUBLIC default)
target_include_directories(test_creeper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
