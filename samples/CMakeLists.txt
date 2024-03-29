file(GLOB samples_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

find_package(OpenMP)

foreach(sample_filename ${samples_list})
  get_filename_component(sample ${sample_filename} NAME_WE)

  add_executable(${sample} ${sample_filename})
  target_link_libraries(${sample} ${PROJECT_LIBRARY} ${PROJECT_LIBRARY_PROBLEMS})
  if(OpenMP_CXX_FOUND)
      target_link_libraries(${sample} OpenMP::OpenMP_CXX)
  endif()
  set_target_properties(${sample} PROPERTIES
    OUTPUT_NAME "${sample}"
    PROJECT_LABEL "${sample}"
    RUNTIME_OUTPUT_DIRECTORY "../"
  )

  # REPORT
  message( STATUS "")
  message( STATUS "General configuration for ${sample}")
  message( STATUS "OK")
  message( STATUS "")

endforeach()