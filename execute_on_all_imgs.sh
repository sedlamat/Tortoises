
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
	#echo $img_dir_out$img_out
	./tortoise $img_dir_in$img_in $img_dir_out$img_out
done
