set(CAPSID_CPP_FILES
    "Capsid.cpp"
    "Utils.cpp"
    "Optimize.cpp"
    "Files.cpp"
    "icosahedron/ico_generator.cpp"
)

set(CAPSID_HPP_FILES 
    "Capsid.h"
    "Utils.h"
    "Optimize.h"
    "Files.h"
    "Defs.h"
    "CSV.h" 
    "Sequence.h"
    "icosahedron/ico_generator.h"
)


include(CheckIPOSupported)

check_ipo_supported(RESULT result OUTPUT output)

add_library(Capsid STATIC ${CAPSID_CPP_FILES} ${CAPSID_HPP_FILES})

if (NOT WIN32)
    target_link_libraries(Capsid PUBLIC m) # -lmath on linux
endif()

# Add source to this project's executable.
add_executable(Phase "main.cpp")
target_link_libraries(Phase Capsid)

if (result)
    set_property(TARGET Capsid PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET Phase PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO is not supported : ${output}")
endif()
