PYTHONPATH=$PYTHONPATH:/mnt/EXT/Schloss-data/kiverson/MG-RAST-Tools/tools/lib/
MG_DOWNLOAD='/mnt/EXT/Schloss-data/kiverson/MG-RAST-Tools/tools/bin/mg-download.py'

# get a list of the metagenomes. The sed bits clean up the text so the end
# result is a list of the metagenomes
$MG_DOWNLOAD --project mgp265 --list | sed 's/\t.*$//g' | sed '1d' | uniq > metagenomes.list


# download the raw reads from mg-rast using their api url. Can be done in
# parallel if you have gnu parallel on your system. If not, uncomment this loop.
#
#for MG in metagenomes.list; do
#	curl http://api.metagenomics.anl.gov//download/${MG}?file=050.1 | tee >(md5sum > {}.md5sum)  | gzip > ${MG}.fna.gz
#done

mkdir data
cat metagenomes.list | parallel "curl http://api.metagenomics.anl.gov//download/{}?file=050.1 | tee >(md5sum > data/{}.md5sum)  | gzip > data/{}.fna.gz"
