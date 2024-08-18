# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


# Only build if we have qhull
if(TARGET external-qhull)

	set(SOURCES filter_close_holes.cpp)

	set(HEADERS filter_close_holes.h)

	add_meshlab_plugin(filter_close_holes2 ${SOURCES} ${HEADERS})

	target_link_libraries(filter_close_holes2 PRIVATE external-qhull)

else()
	message(STATUS "Skipping filter_qhull - missing qhull")
endif()
