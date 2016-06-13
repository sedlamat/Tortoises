
img_dir_in=$HOME"/Images/Tortoises/"
img_dir_out=$img_dir_in"Results/"
img_dir_check=$img_dir_out"GOOD/"

if [ ! -d $img_dir_out ]
then
	mkdir $img_dir_out
fi

imgs=$(ls -r $img_dir_in | grep ".jpg")
#echo $imgs
for img_in in $imgs
do
	echo "Processing image $img_in"
	img_out=res${img_in:0:7}".jpg"

	img_in_path=$img_dir_in$img_in
	img_out_path=$img_dir_out$img_out
	img_check_path=$img_dir_check$img_out

	if [ ! -f $img_check_path ]
	then
		./tortoise $img_in_path $img_out_path
	fi

done
