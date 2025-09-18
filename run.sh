numtask=5
folder_path2="Result/tmp"
folder_path="Result"
mesh_folder_path="Result/tmp/mesh"
fwd_folder_path="Result/tmp/fwd"
exe="H2DMT"
#delete tmp
if [ -d "$folder_path" ]; then
    find "$folder_path" -mindepth 1 -delete
    echo "Folder $folder_path cleaned."
else
    echo "Folder $folder_path does not exist."
fi

mkdir $folder_path2
mkdir $mesh_folder_path
mkdir $fwd_folder_path

mpirun -np $numtask ./$exe
#./$exe
