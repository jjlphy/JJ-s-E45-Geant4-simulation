#!/bin/sh

set -e

main_dir=$(dirname `readlink -f $0`)
src_dir=$main_dir/src
build_dir=$main_dir/.build
bin_dir=$main_dir/bin

rm -rf $build_dir $bin_dir $src_dir/Dict.cc $src_dir/Dict_rdict.pcm
