
img_dir_in=$HOME"/Images/Tortoises/"
img_dir_out=$img_dir_in"Results/"

if [ ! -f $img_dir_out ]
then
	mkdir $img_dir_out
fi

imgs=$(ls $img_dir_in | grep "*.pnm")

for img_in in $imgs
do
	echo "Processing image $img_in"
	img_out="res"${img_in:0:7}".jpg"
	#echo $img_dir_out$img_out
	./tortoise $img_dir_in$img_in $img_dir_out$img_out
done
