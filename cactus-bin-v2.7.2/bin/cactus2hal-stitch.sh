#!/bin/sh
# Hacky shell script to automate stitching two alignments together
set -o errexit
set -o nounset
set -o pipefail
MAX_PARALLELISM=20

if [ $# -ne 5 ]; then
    echo "Stitches together a supertree and subtree alignment. The supertree must have the stitch point as a leaf."
    echo "Usage: $0 supertree_hal subtree_work_dir stitchpoint_name_in_supertree stitchpoint_name_in_subtree destination_hal"
    exit 1
fi
SUPERTREE_HAL=$(readlink -f "$1")
SUBTREE_WORKDIR=$(readlink -f "$2")
STITCHPOINT_NAME1=$3
STITCHPOINT_NAME2=$4
DEST_HAL=$5

echo "Copying supertree to $DEST_HAL"
cp "$SUPERTREE_HAL" "$DEST_HAL"
STITCH_WORKDIR=$(mktemp -d)
echo "Copying subtree workdir to $STITCH_WORKDIR"
trap "rm -fr $STITCH_WORKDIR; exit 1" EXIT
cp -r "$SUBTREE_WORKDIR"/* "$STITCH_WORKDIR"

echo "Substituting names"
# change progressiveAlignment project XML
sed -ir "s/name=\"$STITCHPOINT_NAME2\"/name=\"$STITCHPOINT_NAME1\"/g;s^$SUBTREE_WORKDIR^$STITCH_WORKDIR^g;s/$STITCHPOINT_NAME2\([:,;]\)/$STITCHPOINT_NAME1\1/g" "$STITCH_WORKDIR"/progressiveAlignment/progressiveAlignment_project.xml
# change experiment XMLs
find "$STITCH_WORKDIR" -type f | grep -E '_experiment.xml$' | xargs -P $MAX_PARALLELISM -n 1 sed -i "s/$STITCHPOINT_NAME2;/$STITCHPOINT_NAME1;/g;s^$SUBTREE_WORKDIR^$STITCH_WORKDIR^g"
# change c2h files
find "$STITCH_WORKDIR" -type f | grep -E '.c2h$' | xargs -P $MAX_PARALLELISM -n 1 sed -i "s/s	'$STITCHPOINT_NAME2'/s	'$STITCHPOINT_NAME1'/g"
echo "Appending subtree"
cactus2hal.py --append --event $STITCHPOINT_NAME1 "$STITCH_WORKDIR"/progressiveAlignment/progressiveAlignment_project.xml "$DEST_HAL"

rm -fr "$STITCH_WORKDIR"
