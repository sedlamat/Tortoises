
img_dir_in=$HOME"/Images/Tortoises/"
img_dir_out=$img_dir_in"Results/"

if [ ! -d $img_dir_out ]
then
	mkdir $img_dir_out
fi

imgs=$(ls -r $img_dir_in | grep ".jpg")
echo $imgs
for img_in in $imgs
do
	echo "Processing image $img_in"
	img_out=res${img_in:0:7}".jpg"

	img_in=$img_dir_in$img_in
	img_out=$img_dir_out$img_out
	#echo $img_dir_out$img_out

	if [ ! -f $img_out ]
	then
		./tortoise $img_in $img_out
	fi
	
done
