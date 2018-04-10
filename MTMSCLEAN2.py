# Deconvolution using CLEAN algorithm
import glob
import os


prefix = '3c84-north'
freqs = [1311, 1330, 1350, 1370, 1392, 1410, 1428, 1450]
#measurement sets in list
L = glob.glob('3c84-north.*.MS')
L.sort()

# assign list FL to phase rotated measurements
FL = ['3c84-north.'+str(fr)+'_fixed.MS' for fr in freqs]

default('clean')
'''
###### I: dirty image #####
dirtyname = prefix + '_dirtyimg'

clean(vis = L, imagename = dirtyname, imagermode = 'csclean', psfmode = 'clark',
niter = 0, mode = 'mfs', weighting='uniform', cell = ['2arcsec','2arcsec'],
imsize = [4096,4096], interactive=False)
#dirty image
dirty_image = dirtyname +'.image'

# Export casa dirty image to fits image
exportfits(imagename = dirty_image,fitsimage = dirty_image+'.fits')
#view dirty image statistics
imstat(dirty_image)

################################################################################
################ II: Deconvolve #################################################

cleanname= prefix+'_multiscale_clean' #_round2_60s' #_30s_acal'
#6.2761097215625
clean(vis = L, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '1mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84-north_clean.mask')
imstat(cleanname+'.model')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)

################ III fixvis around 3c84 #########################################


for i, j in enumerate(L):
    #
    freq = freqs[i]
    print "fixing freq = %s, MS = %s" %(freq, j)
    fixvis(vis=j, outputvis='3c84-north.'+str(freq)+'_fixed.MS', phasecenter='J2000 3h19m48.1s 41d30m42s')

################ IV: Deconvolve fixed MSs #################################################

cleanname= prefix+'_multiscale_fixed'
#6.2761097215625
clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '1mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84-north_clean_fixed.mask')
imstat(cleanname+'.model')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)


################ V: Self-calibrate fixed #################################################

for i, j in enumerate(FL):
    #ft(vis=j, model='3c84-north_clean.model', usescratch=True) #_resid.model', usescratch=True) #image model
    #ft(vis=j, model='only_3c84', usescratch=True) #_resid.model', usescratch=True) #image model
    #Bandpass calibration
    freq = freqs[i]
    band = '/net/tina/vault2-tina/mungwariri/17/3c48.'+str(freq)+'.B0'
    calt = 'caltable_'+str(freq)+'_120s.G0'
    print "calibrating 3c84-north.%s.MS using %s" %(freq, band)
    gaincal(vis=j, gaintype='G', solint='120s',combine='scan', gaintable=band, caltable=calt, calmode='p')
    applycal(vis=j, gaintable=[band, calt])
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna',
    figfile='selfcal_'+str(freq)+'_phaseVtime_120s_pcal.png', showgui=False)
######################################################################################################################


################ VI: Deconvolve selfcaled fixed MSs #################################################
cleanname= prefix+'_multiscale_fixed_120s' #_round2_60s' #_30s_acal'
#10.0417755mJy
clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '1mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84-north_multiscale_fixed.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)



################ VII: Mask the model image (cleanname_model) #####################################
# In IPython notebook

################ VIII: Subtract sources besides 3c84 from field ########################
for j in FL:
    ft(vis=j, model='no3c84_masked.model', usescratch=True)
    uvsub(vis = j)

# os.system('cd only_3c84/')

################ IX: Deconvolve, 3c84 only ########################

cleanname= '3c84_only_subbed' #_round2_60s' #_30s_acal'
#10.0417755mJy
clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '1mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)



###########################################################################################

################ X: Selfcal 2, 3c84 only ########################
for i, j in enumerate(FL):
    mod = '3c84_masked.model' #model image containing 3c84 (other sources blanked)
    ft(vis=j, model=mod, usescratch=True) #insert this model into the MS
    print 'inserted %s into %s' %(mod, j)
    ####
    #replace the DATA column with the CORRECTED_DATA column
    tb.open(j,nomodify=False)                  #open MS
    corrected = tb.getcol('CORRECTED_DATA')    # read the CORRECTED_DATA column #i.e calibrated visibilities
    data = tb.getcol('DATA')                   # read the DATA column
    tb.putcol('DATA',corrected)                # now the DATA column contains the calibrated visibilities
    print 'inserted CORRECTED_DATA column into DATA column for', j
    ###
    #self calibrate
    freq = freqs[i]

    calt = 'caltable_'+str(freq)+'_300s.G0'       #calibration table
    print "self calibrating 3c84-north.%s.MS" %freq
    gaincal(vis=j, gaintype='G', solint='300s',combine='scan', caltable=calt, calmode='p') #phase-only gain calibration
    print "applying gain calibration only!"
    applycal(vis=j, gaintable= calt)   #apply gain calibration  | plot result
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna', figfile='selfcal_'+str(freq)+'_phaseVtime_300s_pcal.png', showgui=False)

    tb.putcol('DATA',data)                      # now the DATA column contains the raw data again
    tb.close()
    tb.done()                           #Done restoring MS vis columns
    print 'done with', j        #move on to next MS or exit

###########################################################################################
################ XI: Deconvolve, 3c84 only ########################

cleanname= '3c84_only_subbed_300s' #_round2_60s' #_30s_acal'

clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '1mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)

###########################################################################################
'''

# os.system('cd second_120s')

################ XII: Selfcal3, shorter solint = 120s, 3c84 only ########################
for i, j in enumerate(FL):
    mod = '3c84_only_subbed_300s.model.tt0' #model image containing 3c84 (other sources blanked)
    ft(vis=j, model=mod, usescratch=True)   #insert this model into the MS
    print 'inserted %s into %s' %(mod, j)
    ####
    #replace the DATA column with the CORRECTED_DATA column
    tb.open(j,nomodify=False)                  #open MS
    corrected = tb.getcol('CORRECTED_DATA')    # read the CORRECTED_DATA column #i.e calibrated visibilities
    data = tb.getcol('DATA')                   # read the DATA column
    tb.putcol('DATA',corrected)                # now the DATA column contains the calibrated visibilities
    print 'inserted CORRECTED_DATA column into DATA column for', j
    ###
    #self calibrate
    freq = freqs[i]

    calt = 'caltable_'+str(freq)+'_120s.G0'       #calibration table
    print "self calibrating 3c84-north.%s.MS" %freq
    gaincal(vis=j, gaintype='G', solint='120s',combine='scan', caltable=calt, calmode='p') #phase-only gain calibration
    print "applying gain calibration only!"
    applycal(vis=j, gaintable= calt)   #apply gain calibration  | plot result
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna', figfile='selfcal_'+str(freq)+'_phaseVtime_120s_pcal.png', showgui=False)

    tb.putcol('DATA',data)        # now the DATA column contains the raw data again
    tb.close()
    tb.done()                     #Done restoring MS vis columns
    print 'done with', j          #move on to next MS or exit

###########################################################################################
'''################ XIII: Deconvolve, 3c84 only thr = 0.9mJy ########################'''

cleanname= '3c84_only_subbed_120s'

clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '0.9mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)


###########################################################################################
'''################ XIV: Selfcal4, shorter solint = 60s, 3c84 only ########################'''
for i, j in enumerate(FL):
    mod = '3c84_only_subbed_120s.model.tt0' #model image containing 3c84 (other sources blanked)
    ft(vis=j, model=mod, usescratch=True) #insert this model into the MS
    print 'inserted %s into %s' %(mod, j)
    ####
    #replace the DATA column with the CORRECTED_DATA column
    tb.open(j,nomodify=False)                  #open MS
    corrected = tb.getcol('CORRECTED_DATA')    # read the CORRECTED_DATA column #i.e calibrated visibilities
    data = tb.getcol('DATA')                   # read the DATA column
    tb.putcol('DATA',corrected)                # now the DATA column contains the calibrated visibilities
    print 'inserted CORRECTED_DATA column into DATA column for', j
    ###
    #self calibrate
    freq = freqs[i]

    calt = 'caltable_'+str(freq)+'_60s.G0'       #calibration table
    print "self calibrating 3c84-north.%s.MS" %freq
    gaincal(vis=j, gaintype='G', solint='60s',combine='scan', caltable=calt, calmode='p') #phase-only gain calibration
    print "applying gain calibration only!"
    applycal(vis=j, gaintable= calt)   #apply gain calibration  | plot result
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna', figfile='selfcal_'+str(freq)+'_phaseVtime_60s_pcal.png', showgui=False)

    tb.putcol('DATA',data)        # now the DATA column contains the raw data again
    tb.close()
    tb.done()                     #Done restoring MS vis columns
    print 'done with', j          #move on to next MS or exit

###########################################################################################
'''################ XV: Deconvolve, 3c84 only thre = 0.8mJy ########################'''

cleanname= '3c84_only_subbed_60s' #_round2_60s' #_30s_acal'

clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '0.8mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)



##############################################################################################
'''################ XVI: Selfcal5, shorter solint = 30s, 3c84 only ########################'''
for i, j in enumerate(FL):
    mod = '3c84_only_subbed_60s.model.tt0' #model image containing 3c84 (other sources blanked)
    ft(vis=j, model=mod, usescratch=True) #insert this model into the MS
    print 'inserted %s into %s' %(mod, j)
    ####
    #replace the DATA column with the CORRECTED_DATA column
    tb.open(j,nomodify=False)                  #open MS
    corrected = tb.getcol('CORRECTED_DATA')    # read the CORRECTED_DATA column #i.e calibrated visibilities
    data = tb.getcol('DATA')                   # read the DATA column
    tb.putcol('DATA',corrected)                # now the DATA column contains the calibrated visibilities
    print 'inserted CORRECTED_DATA column into DATA column for', j
    ###
    #self calibrate
    freq = freqs[i]

    calt = 'caltable_'+str(freq)+'_30s_apcal.G0'       #calibration table
    print "self calibrating 3c84-north.%s.MS" %freq
    gaincal(vis=j, gaintype='G', solint='30s',combine='scan', caltable=calt, calmode='p') #amp+phase gain calibration
    print "applying amp and phase gain calibration only!"
    applycal(vis=j, gaintable= calt)   #apply gain calibration  | plot result
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna', figfile='selfcal_'+str(freq)+'_phaseVtime_30s_apcal.png', showgui=False)

    tb.putcol('DATA',data)        # now the DATA column contains the raw data again
    tb.close()
    tb.done()                     #Done restoring MS vis columns
    print 'done with', j          #move on to next MS or exit

##########################################################################################
'''################ XVII: Deconvolve, 3c84 only thre = 0.75mJy ########################'''

cleanname= '3c84_only_subbed_30s' #_30s_pcal'

clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '0.75mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)

#Last step Amplitude cal

##############################################################################################
'''################ XVI: Selfcal6: ampcal, solint = inf, 3c84 only ########################'''
for i, j in enumerate(FL):
    mod = '3c84_only_subbed_30s.model.tt0' #model image containing 3c84 (other sources blanked)
    ft(vis=j, model=mod, usescratch=True) #insert this model into the MS
    print 'inserted %s into %s' %(mod, j)
    ####
    #replace the DATA column with the CORRECTED_DATA column
    tb.open(j,nomodify=False)                  #open MS
    corrected = tb.getcol('CORRECTED_DATA')    # read the CORRECTED_DATA column #i.e calibrated visibilities
    data = tb.getcol('DATA')                   # read the DATA column
    tb.putcol('DATA',corrected)                # now the DATA column contains the calibrated visibilities
    print 'inserted CORRECTED_DATA column into DATA column for', j
    ###
    #self calibrate
    freq = freqs[i]

    calt = 'caltable_'+str(freq)+'_acal.G0'       #calibration table
    print "self calibrating 3c84-north.%s.MS" %freq
    gaincal(vis=j, gaintype='G', solint='inf',combine='scan', caltable=calt, calmode='a') #amp gain calibration
    print "applying amp and phase gain calibration only!"
    applycal(vis=j, gaintable= calt)   #apply gain calibration  | plot result
    plotcal(caltable=calt, xaxis='time', yaxis='phase', poln='X', subplot =221, iteration = 'antenna', figfile='selfcal_'+str(freq)+'_phaseVtime_acal.png', showgui=False)

    tb.putcol('DATA',data)        # now the DATA column contains the raw data again
    tb.close()
    tb.done()                     #Done restoring MS vis columns
    print 'done with', j          #move on to next MS or exit

##########################################################################################
'''################ XVII: Deconvolve, 3c84 only thre = 0.75mJy ########################'''

cleanname= '3c84_only_subbed_ampcal'

clean(vis = FL, imagename = cleanname, imagermode = 'csclean', psfmode = 'clark',
threshold = '0.75mJy', niter = 20000, mode = 'mfs', weighting='uniform',
cell = ['2arcsec','2arcsec'], multiscale=[5,10,13], smallscalebias=0.01, nterms= 2,
imsize = [4096,4096], interactive=False, mask='3c84_only.mask')
imstat(cleanname+'.model.tt0')

cleanname_restored = cleanname+ '.image.tt0'
# exportfits
exportfits(imagename = cleanname_restored,fitsimage = cleanname_restored+'.fits', overwrite = True)

# model image:
cleanname_model = cleanname+'.model.tt0'
# exportfits
exportfits(imagename = cleanname_model,fitsimage = cleanname_model+'.fits', overwrite = True)

# residual image:
cleanname_residual = cleanname+'.residual.tt0'
# exportfits
exportfits(imagename = cleanname_residual,fitsimage = cleanname_residual+'.fits', overwrite = True)

# synthesized (dirty) beam
cleanname_psf = cleanname+'.psf.tt0'
# exportfits
exportfits(imagename = cleanname_psf,fitsimage = cleanname_psf+'.fits', overwrite = True)
