#! /bin/bash
#Create directories and busco fastas with clean IDs 
d=$1
nd=$(echo $d | sed 's/\./_/g' | sed 's/.busco//g' | sed 's/-/_/g') 
echo $nd
mkdir $nd 
cd $nd
cp ../$d/single_copy_busco_sequences/* .  
sed -i "s/^>.*/>$nd/g" *
cd ..

