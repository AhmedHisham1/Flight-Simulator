#! /bin/csh 

cd /usr/share/games/flightgear

#setenv LD_LIBRARY_PATH /usr/share/games/flightgear/lib:$LD_LIBRARY_PATH
#setenv FG_ROOT /usr/share/games/flightgear/data
#setenv FG_SCENERY /usr/share/games/flightgear/Scenery:$FG_ROOT/Scenery:$FG_ROOT/WorldScenery

fgfs --aircraft=dhc2F --fdm=network,localhost,5501,5502,5503 --fog-fastest --disable-clouds --start-date-lat=2004:06:01:09:00:00 --disable-sound --in-air --enable-freeze --airport=KSFO --runway=10L --altitude=7224 --heading=113 --offset-distance=4.72 --offset-azimuth=0
