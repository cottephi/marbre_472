#! /bin/bash

if [ $# -ne 1 ] ; then
   echo "ERROR: need one argument, the folder to merge"
   exit 1
fi

path="$1"
if [ ! -d $path ] ; then
   echo "$path : No such directory"
   exit 1
fi

echo "Sorting files in $path..."

if ! ls $path/*r*_h* 2>/dev/null 1>&2 ; then
  echo "No data to sort (will not sort if only marble files are present)"
  exit 1
fi

if ! ls $path/*marble* 2>/dev/null 1>&2 ; then
  echo "No marble to sort (will not sort if only data files are present)"
  exit 1
fi

if [ -d $path/data ] ; then
  echo "Removing previous data files..."
  rm $path/data/*
else
  mkdir $path/data
fi
if [ -d $path/cali ] ; then
  echo "Removing previous cali files..."
  rm $path/cali/*
else
  mkdir $path/cali
fi

echo "Sorting data..."
mv $path/*marble* $path/cali/
mv $path/*r*_h* $path/data/

sleep 1

if  [[ $(find ${path} -maxdepth 0 -type f | wc -l) -gt 0 ]] ; then 
  if [ -d $path/other_files ] ; then
    echo "Removing previous other files..."
    rm $path/other_files/*
  else
  mkdir $path/other_files
  fi
  for File in $path/*; do
   if ! [ -d $File ]; then
     mv $File $path/other_files/
   fi
  done
fi

echo "...done"

./merge.sh $path/data
./merge.sh $path/cali
