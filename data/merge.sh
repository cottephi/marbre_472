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
echo "merging file in $path..."

if ls $path/*merged* 2>/dev/null 1>&2 ; then
  rm $path/*merged*
  echo "removed previous merged file..."
fi

for file in $path/*_*1_*.csv ; do
  filename="$(basename "${file}")"
  if [ ! -f ${file} ] ; then
    continue
  fi
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/merged_1.csv
done
for file in $path/*_*2_*.csv ; do
  filename="$(basename "${file}")"
  if [ ! -f ${file} ] ; then
    continue
  fi
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/merged_2.csv
done
for file in $path/*_*3_*.csv ; do
  filename="$(basename "${file}")"
  if [ ! -f ${file} ] ; then
    continue
  fi
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/merged_3.csv
done
for file in $path/*_*4_*.csv ; do
  filename="$(basename "${file}")"
  if [ ! -f ${file} ] ; then
    continue
  fi
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/merged_4.csv
done
for file in $path/*_*5_*.csv ; do
  filename="$(basename "${file}")"
  if [ ! -f ${file} ] ; then
    continue
  fi
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/merged_5.csv
done
echo "done"
exit 0
