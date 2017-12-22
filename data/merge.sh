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

for file in $path/*_r1_*1000Hz.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r1_1000Hz_merged.csv
done
for file in $path/*_r2_*1000Hz.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r2_1000Hz_merged.csv
done
for file in $path/*_r3_*1000Hz.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r3_1000Hz_merged.csv
done
for file in $path/*_r4_*1000Hz.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r4_1000Hz_merged.csv
done
for file in $path/*_r5_*1000Hz.csv ; do
  filename="$(basename "${file}")"
  echo "Merging $filename..."
  tail -n +2 "$file" | cat >> $path/r5_1000Hz_merged.csv
done
echo "done"
exit 0
