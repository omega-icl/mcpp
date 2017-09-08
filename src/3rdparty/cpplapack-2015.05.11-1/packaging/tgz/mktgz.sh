#!/bin/sh

version=2015.05.11-1

name=cpplapack-$version
dir=/tmp/$name

mkdir $dir || exit 1
cp -r ../../* $dir || exit 1
rm -rf `find $dir -type d -name .svn` || exit 1

cd /tmp || exit 1
rm -f $name.tar.gz || exit 1
tar czf $name.tar.gz $name || exit 1
rm -rf $dir || exit 1
mv $name.tar.gz ~/ || exit 1
echo "~/$name.tar.gz was successfully created."
