#!/usr/bin/bash

#
# automation script to package both Debug & Release builds with CPack
#

# note: could be extended for cross-compiling

# tested on:
#   - MinGW64 (Windows) : MSVC 2022
#   - WSL (Ubuntu) : gcc 11 / clang 19

project_dir="${PWD}"
build_dir="_build"
inst_dir="_inst"

cmake_cmd=`command -v "cmake"`
cpack_cmd=`command -v "cpack"`

message () {
    local msg=$1
    local green_color="\e[92m"
    local end_color="\e[0m"
    echo -e "${green_color}== ${msg} ==${end_color}"
}

cmake_get_generator () {
    if [[  "$OSTYPE" == "win32" || "$OSTYPE" == "msys" ]]; then
        local cmake_generator="'Visual Studio 17 2022' -A x64"
    elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "cygwin" || "$OSTYPE" == "freebsd"* ]]; then
        local cmake_generator="'Unix Makefiles'"
    fi
    echo "${cmake_generator}"
}

cmake_get_build_options () {
    if [[  "$OSTYPE" == "win32" || "$OSTYPE" == "msys" ]]; then
        local cmake_build_options="-- -v:minimal -m:16"
    elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "cygwin" || "$OSTYPE" == "freebsd"* ]]; then
        local cmake_build_options="-- -j 16"
    fi
    echo "${cmake_build_options}"
}

cmake_configure_and_build () {

    local build_type=$1

    message "Configure: ${build_type}"

    local cmake_generator="$(cmake_get_generator)"

    local cmd="\"$cmake_cmd\" -G ${cmake_generator} "
    cmd+="-D \"CMAKE_INSTALL_PREFIX:string=../../${inst_dir}\" "
    if [[ "$OSTYPE" == "win32" || "$OSTYPE" == "msys" ]]; then
        cmd+="-D \"CMAKE_CONFIGURATION_TYPES:string=${build_type}\" "
    else
        cmd+="-D \"CMAKE_BUILD_TYPE:string=${build_type}\" "
    fi
    cmd+="../.."
    echo $cmd
    eval $cmd

    message "Build: ${build_type}"

    local cmake_build_options="$(cmake_get_build_options)"

    local cmd="\"$cmake_cmd\" --build . --config ${build_type} --target install ${cmake_build_options}"
    echo $cmd
    eval $cmd
}

cmake_generate_package_config_file () {

    local build_types=("$@")

    cpack_multi_config_file="${project_dir}/${build_dir}/CPackMultiConfig.cmake"

    if [ ! -f "${cpack_multi_config_file}" ]; then

        echo "include(\"${build_types[-1]}/CPackConfig.cmake\")" > ${cpack_multi_config_file}
        echo "set(CPACK_INSTALL_CMAKE_PROJECTS" >> ${cpack_multi_config_file}            
        for build_type in "${build_types[@]}"
        do
            echo "    \"${build_type};OSD_lite;All;/\"" >> ${cpack_multi_config_file}
        done
        echo ")" >> ${cpack_multi_config_file}
    fi

    echo ${cpack_multi_config_file}
}

cmake_package () {

    local build_types=("$@")

    message "Package"

    cd "${project_dir}/${build_dir}/"

    local cpack_multi_config_file="$(cmake_generate_package_config_file ${build_types[@]})"

    local cmd="\"$cpack_cmd\" --config \"${cpack_multi_config_file}\""
    echo $cmd
    eval $cmd

    cd "${project_dir}"
}

main () {

    if [ ! -d "${project_dir}/${build_dir}" ]; then
        mkdir "${project_dir}/${build_dir}"
    fi

    local build_types=("Debug" "Release")
    for build_type in "${build_types[@]}"
    do
        if [ ! -d "${project_dir}/${build_dir}/${build_type}" ]; then
            mkdir "${project_dir}/${build_dir}/${build_type}"
        fi

        cd "${project_dir}/${build_dir}/${build_type}"

        cmake_configure_and_build "${build_type}"
    done

    cmake_package "${build_types[@]}"
}

main "$@"

