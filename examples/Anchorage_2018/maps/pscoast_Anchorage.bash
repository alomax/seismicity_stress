
OUTFILE=Anchorage_2018.xy

# seismicity.grid.latlontext=38.55,38.8,-28.35,-28.0,0.4,4.0

# GMT 4
cat << END > ${OUTFILE}
>       GMT_LONLAT                                                                              
END
pscoast -Jx1d -R-152.4/-147.5/60.1/62.7 -W -Df -Ia -M >> ${OUTFILE}



# GMT 6
cat << END > ${OUTFILE}
>       GMT_LONLAT                                                                              
END
gmt pscoast -Jx1d -R-152.4/-147.5/60.1/62.7 -W -Df -M >> ${OUTFILE}
gmt pscoast -Jx1d -R-152.4/-147.5/60.1/62.7 -Df -Ia -M >> ${OUTFILE}



