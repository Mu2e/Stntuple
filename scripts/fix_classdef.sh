#1/usr/bin/bash

flist=`find  Stntuple -name \*.hh | xargs grep ClassDef | grep -v ClassDefOverride | awk -F : '{print $1}'`
for f in $flist ; do
    echo $f
    cat $f | sed 's/ClassDef/ClassDefOverride/' >| $f.1
    mv $f.1 $f
done
     
