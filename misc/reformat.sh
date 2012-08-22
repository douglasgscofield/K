#!/bin/sh

for f in *.cpp *.h ; do
    origf="${f}.orig";
    newf="${f}.new";
    cp -f $f $origf
    cat $f | sed -e 's/\/\*\/\//\/\* \//g' -e 's/\/\/\*\//\/ \*\//g' -e 's/\r$//g' | indent -cdw -bls -nbad -bap -nbbo -nbc -br -c33 -cd33 -ncdb -ce -ci4 -cli0 -cp33 -cs -d0 -di1 -nfc1 -nfca -nhnl -i4 -ip0 -l79 -lp -npcs -nprs -psl -saf -sai -saw -nsc -nsob -nss -ut -ts4 -bfda > $newf
    mv -f $newf $f
done
