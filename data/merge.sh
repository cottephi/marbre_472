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

for file in $path/*_r1_*.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r1_merged.csv
done
for file in $path/*_r2_*.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r2_merged.csv
done
for file in $path/*_r3_*.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r3_merged.csv
done
for file in $path/*_r4_*.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r4_merged.csv
done
for file in $path/*_r5_*.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r5_merged.csv
done
echo "done"
exit 0
