@PACKAGE_INIT@

set("@PROJECT_NAME@_WITH_CBOOL" @WITH_CBOOL@)
set("@PROJECT_NAME@_WITH_QP" @WITH_QP@)
set("@PROJECT_NAME@_WITH_XDP" @WITH_XDP@)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()
