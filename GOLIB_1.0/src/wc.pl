#!/usr/bin/perl -w
# line count of the MC++ project

$src1Directory="mc";
$src2Directory="global";

$directory[0]="$src1Directory/*.h $src1Directory/*.cpp";
$directory[1]="$src2Directory/*.h $src2Directory/*.cpp";

system "wc -l @directory";
