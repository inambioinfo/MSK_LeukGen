#!/bin/bash  
# usage bash get_module_version.sh MODULENAME
# module name as first argument
mod=$1
# get version
perl -M$mod -e "print \"$mod: \$$mod::VERSION\n\""