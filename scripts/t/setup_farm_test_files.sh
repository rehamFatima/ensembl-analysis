 #!/bin/bash

farm_test_dir=$EAT_DIR'/ensembl-analysis/modules/t/farm'

#read in values
db_user=$1  # for EAT_DB_USER
db_pass=$2  # for EAT_DB_PASS
db_name=$3  # for EAT_DB_NAME
output_path=$4  # for EAT_OUTPUT_PATH
db_user_ro=$5  # for EAT_DB_USER_RO
db_host=$6  # for EAT_DB_HOST

#get a list of template files and update
echo Working in $farm_test_dir
for f in  $farm_test_dir/*.template ;
do
 filename=${f##*/}
 output_file=${filename%.*}
 echo Using template $filename to generate $output_file
 gawk '{ gsub("EAT_DB_USER", "$db_user"); print }' | \
     gawk '{ gsub("EAT_DB_PASS", "$db_pass"); print }' | \
     gawk '{ gsub("EAT_DB_NAME", "$db_name"); print }' | \
     gawk '{ gsub("EAT_DB_USER", "$db_user_ro"); print }' | \
     gawk '{ gsub("EAT_DB_HOST", "$db_host"); print }' > $farm_test_dir/$output_file
done

echo Working in $farm_test_dir/conf
for f in  $farm_test_dir/conf/*.template ;
do
 filename=${f##*/}
 output_file=${filename%.*}
 echo Using template $filename to generate $output_file
done
