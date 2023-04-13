for dir in `find . -type d`; do [ "$dir" != "." ] && mkdir -p $BASILISK/skele/$dir ;done
for file in `find . -type f`; do ln -s $SKELE_SRC/$file $BASILISK/skele/$file; done
