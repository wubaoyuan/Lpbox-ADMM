echo "Starting testing on glass task"
./clustering --input clustering_data/glass.txt --sample 214 --class 7 --dim 9 --neighbor 5 --permute 0 --log glass_out.txt;
echo "Glass task finishes"
echo
echo "Starting testing on iris task"
./clustering --input clustering_data/iris.txt --sample 150 --class 3 --dim 4 --neighbor 5 --permute 1 --log iris_out.txt;
echo "Iris task finishes"
echo
echo "Starting testing on wine task"
./clustering --input clustering_data/wine.txt --sample 178 --class 3 --dim 13 --neighbor 5 --permute 2 --log wine_out.txt;
echo "Wine task finishes"
