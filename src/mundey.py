import sys
import numpy as np 
from astropy.io import fits
import xml.etree.ElementTree as ET
import warnings
import lightkurve
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.units as u
from IPython.core.display import display, HTML

from tqdm import tqdm

import pkg_resources


class mundey_tpf(lightkurve.TessTargetPixelFile):


# =========================================================================
# Main routine to perform a new calibration of a TESS target pixel file
# =========================================================================

	def calibrate(self,verbose=True,smear='automatic',smear_range=None,wings=None,ddir='./',target=None,use_simbad=True):
		"""Performs a new calibration for a target pixel file
			Parameters
			----------
			smear : 'automatric', standard', or 'alternate'
				If the string 'automatic' is passed, will automatically work out which smear correction to use
				If the string 'standard' is passed, will calculate the smear correction directly from the collateral files.
				If the string 'alternate' is passed, smear correction will be estimated from the target pixel file.
				If the string 'both' is passed, smear correction will be estimated from both the and the target pixel file.

			Returns
			-------

		"""

		#Identify star we're looking at
		if target is not None:
			target = target 
		elif use_simbad:
			customSimbad = Simbad()
			customSimbad.add_votable_fields('flux(V)')

			result_table = customSimbad.query_region(coord.SkyCoord(self.ra,self.dec, unit=(u.deg, u.deg), frame='icrs'))
			result_table.sort(['FLUX_V'])
			target = result_table[0]["MAIN_ID"]
		else:
			target = 'TIC ' + str(self.targetid)
		self.ddir = ddir


		if verbose:
			print("Beginning calibration of "+target)

		#Set up some arrays 
		flux = np.copy(self.hdu[1].data['RAW_CNTS'])
		cols = np.arange(self.column,self.column+self.shape[2])

		#Find the output for each column
		outputs=[]
		for col in cols:
			if (col >=45) & (col <=556):
				outputs.append('A')
			elif (col >= 557) & (col <= 1068):
				outputs.append('B')
			elif (col >= 1069) & (col <= 1580):
				outputs.append('C')
			elif (col >= 1581) & (col <= 2092):
				outputs.append('D')

		path = self.path

		#Define x,y and dimensions of target pixel file
		x1 = self.hdu[2].header['CRVAL1P']-1
		dx = self.hdu[2].header['NAXIS1']

		y1 = self.hdu[2].header['CRVAL2P']-1
		dy = self.hdu[2].header['NAXIS2']

		#Determine which smear correction to use if set to auto
		if smear == 'automatic':
			if y1+dy >= 2048:
				smear = 'alternate'
			else:
				smear = 'standard'

		if verbose:
			print("Sector: "+str(self.sector))
			print("Camera: " +str(self.camera))
			print("CCD: " +str(self.ccd))
			print("Outputs: " +str(set(outputs)))
			print("")
			print("Searching for collateral target pixel files")

		#Search for collateral data

		sm_cal = np.zeros((len(set(outputs)),flux.shape[0],10,512))
		tv_cal = np.zeros((len(set(outputs)),flux.shape[0],2078,11))
		calcsmear = np.zeros(len(set(outputs)),dtype=bool)
		calcblack = np.zeros(len(set(outputs)),dtype=bool)


		for idx,output in enumerate(set(outputs)):

			smearfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-smrow-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col-blf.fits'

			try:
				if verbose:
					print("... loading smear file "+smearfile)
				smrow = fits.open(smearfile)
				sm_cal[idx,:,0,:] = smrow[1].data['SMEAR']
				calcsmear[idx]=False
			
			except FileNotFoundError:
				smrowfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-smrow-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col.fits'
			
				if verbose:
					print("... file not found: "+smearfile)
					print("... smear correction will be calculated from raw counts")
					print("... loading smear row file "+smrowfile)
					print("")

				smrow = fits.open(smrowfile)
				sm_cal[idx] = smrow[1].data['SMROW_RAW']
				calcsmear[idx]=True

	
			onedblackfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-tvcol-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col-blf.fits'
			
			try:
				if verbose:
					print("... loading 1D black file "+onedblackfile)
				tvcol = fits.open(onedblackfile)
				tv_cal[idx,:,:,0] = tvcol[1].data['1DBLACK']
				calcblack[idx]=False

			except FileNotFoundError:
				tvcolfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-tvcol-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col.fits'
				if verbose:
					print("... file not found: "+onedblackfile)
					print("... 1D black correction will be calculated from raw counts")
					print("... loading trailing virtual column file "+tvcolfile)
					print("")

				tvcol = fits.open(tvcolfile)
				tv_cal[idx] = tvcol[1].data['TVCOL_RAW']
				calcblack[idx]=True


		# Subtract fixed offset, add mean black
		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux[:,:,indices] = self.fixedoffset_meanblack(flux[:,:,indices],output,verbose=verbose)
			if calcsmear[idx]:
				sm_cal[idx] = self.fixedoffset_meanblack(sm_cal[idx],output,verbose=False)
			if calcblack[idx]:
				tv_cal[idx] = self.fixedoffset_meanblack(tv_cal[idx],output,verbose=False)
	

		# Remove 2D black
		black2dfile = self.ddir+'tess2018323-'+str(self.camera)+'-'+str(self.ccd)+'-2dblack.fits'

		if verbose:
			print("")
			print("Loading 2D black model")
			print("... "+black2dfile)

		
		black2dfits = fits.open(black2dfile)
		black2dimg = black2dfits[1].data

		flux = self.twodblack(flux,black2dimg[y1:y1+dy,x1:x1+dx],verbose=verbose)

		for idx,output in enumerate(set(outputs)):
			if output == 'A':
				sx1 = 44
				sx2 = 556
				tx1 = 2092
				tx2 = 2103
			if output == 'B':
				sx1 = 556
				sx2 = 1068
				tx1 = 2103
				tx2 = 2114
			if output == 'C':
				sx1 = 1068
				sx2 = 1580
				tx1 = 2114
				tx2 = 2125
			if output == 'D':
				sx1 = 1580
				sx2 = 2092
				tx1 = 2125
				tx2 = 2136
			if calcsmear[idx]:
				sm_cal[idx] = self.twodblack(sm_cal[idx],black2dimg[2058:2068,sx1:sx2],verbose=False)
			if calcblack[idx]:
				tv_cal[idx] = self.twodblack(tv_cal[idx],black2dimg[:,tx1:tx2],verbose=False)


		# Calculate 1D black
		for idx,output in enumerate(set(outputs)):
			if calcblack[idx]:
				tv_cal[idx] = self.calconedblack(tv_cal[idx],output,verbose=verbose)


		# Remove 1D black		
		tv_cal = tv_cal[:,:,:,0] #By now only this slice has the calibrated value

		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux = self.onedblack(flux[:,:,indices],tv_cal[idx,:,y1:y1+dy],verbose=verbose)
			if calcsmear:
				sm_cal[idx] = self.onedblack(sm_cal[idx],tv_cal[idx,:,2058:2068],verbose=verbose)


		# Correct for non-linearity and gain
		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]
		
			flux[:,:,indices] = self.linearity_gain(flux[:,:,indices],output,verbose=verbose)
			if calcsmear[idx]:
				sm_cal[idx] = self.linearity_gain(sm_cal[idx],output,verbose=False)


		# Correct for LDE undershoot
		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux[:,:,indices] = self.undershoot(flux[:,:,indices],output,verbose=verbose)
			if calcsmear[idx]:
				sm_cal[idx] = self.undershoot(sm_cal[idx],output,verbose=False)


		# Calculate smear correction
		for idx,output in enumerate(set(outputs)):
			if calcsmear[idx]:
				sm_cal[idx] = self.calcsmear(flux,sm_cal[idx],output,smear=smear,smear_range=smear_range,wings=wings, verbose=verbose)

		# Remove smear correction
		sm_cal = sm_cal[:,:,0,:] #By now only this slice has the calibrated value
		flux = self.smear(flux,sm_cal,output[0],verbose=verbose)

		# Flatfield correction
		flux = self.flatfield(flux,verbose=verbose)

		# Convert to electrons per second
		flux = self.to_electrons_per_second(flux,verbose=verbose)

		self.hdu[1].data['FLUX'] = flux

		print("")
		print("Calibration Complete")

# =========================================================================
# Subtract fixed offset; add mean black
# =========================================================================

	def fixedoffset_meanblack(self,img,output,verbose=True):
		if verbose:
			print("")
			print("Adjusting for fixed offset and mean black - Output " + output)

		img = img - self.hdu[1].header['FXDOFF']
		img = img + self.hdu[1].header['MEANBLC'+output]

		return img

# =========================================================================
# Remove 2D black
# =========================================================================

	def twodblack(self,img,twodblackimg,verbose=True):
		if verbose:
			print("Performing 2-D black correction")

		img = img - twodblackimg

		return img

# =========================================================================
# Calculate 1D black
# =========================================================================

	def calconedblack(self,tvcol,output,verbose=True):
		if verbose:
			print("")
			print("Calculating 1-D black correction")
			print("... this step may take a while;")
			print("... if you haven't already done so")
			print("... you might like to take this opportunity")
			print("... to join your union.")		
			display(HTML('<a href="https://www.australianunions.org.au/join" target="_blank" rel="noreferrer noopener">Join an Australian Union</a> or <a href="https://www.ituc-csi.org/?page=abook" target="_blank" rel="noreferrer noopener">find an International Trade Union Confederation affiliate in your country</a>.'))

		warnings.simplefilter('ignore',np.RankWarning)

		y1 = self.hdu[2].header['CRVAL2P']-1
		dy = self.hdu[2].header['NAXIS2']

		#For each frame, fit a polynomial to the TV col values, choosing the order using the AIC
		nrows = tvcol.shape[1]
		x = np.arange(nrows)/nrows-0.5
		ks = np.arange(20,50)

		black1d = np.zeros((tvcol.shape[0],tvcol.shape[1]))

		for frame in tqdm(np.arange(tvcol.shape[0]),desc='Frame'):
			aics = []
			y = np.nanmean(tvcol[frame,:,:-1],axis=1) #Average over all but the last column, which seems to be off

			for k in ks:
				p = np.polynomial.Polynomial.fit(x,y,k)
				aic = 2*k + nrows*np.log(np.sum((y-p(x))**2))
				aics.append(aic)

			order = ks[np.argmin(np.asarray(aics))]
			p = np.polynomial.Polynomial.fit(x,y,order)

			black1d[frame,:] = p(x)

		#Save calibration data for later use
		col1 = fits.Column(name=self.hdu[1].header['TTYPE1'], format=self.hdu[1].header['TFORM1'], unit=self.hdu[1].header['TUNIT1'], disp=self.hdu[1].header['TDISP1'], array=self.hdu[1].data['TIME'])
		col2 = fits.Column(name=self.hdu[1].header['TTYPE2'], format=self.hdu[1].header['TFORM2'], disp=self.hdu[1].header['TDISP2'], array=self.hdu[1].data['CADENCENO'])
		col3 = fits.Column(name='1DBLACK', format=str(black1d.shape[1])+'D', unit='count', disp='F14.7', dim='('+str(black1d.shape[1])+')', array=black1d)

		cols = fits.ColDefs([col1,col2,col3])
		blackhdu = fits.BinTableHDU.from_columns(cols)

		path=self.path
		onedblackfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-tvcol-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col-blf.fits'
		blackhdu.writeto(onedblackfile)

		tvcol[:,:,0] = black1d

		return tvcol

# =========================================================================
# Remove 1D black
# =========================================================================	

	def onedblack(self,img,onedblackimg,verbose=True):
		if verbose:
			print("Performing 1-D black correction")

		img = img - np.expand_dims(onedblackimg,axis=2)

		return img		

# =========================================================================
# Correct for non-linearity and gain
# =========================================================================

	def linearity_gain(self,img,output,verbose=True):
		if verbose:
			print("")
			print("Correcting for non-linearity and gain")

		# Convert to mean signal per exposure
		factor = 1.

		# Check to see if cosmic ray mitigation is enabled
		if self.hdu[0].header['CRMITEN'] == True:
			# Then two of every N exposures is discarded
			factor = (self.hdu[0].header['CRBLKSZ']-2)/self.hdu[0].header['CRBLKSZ']

		img = img/(factor*self.hdu[1].header['NUM_FRM'])

		# Load linearity model
		lintree = ET.parse(self.ddir+'tess2018143203310-41006_100-linearity.xml')
		linroot = lintree.getroot()

		poly = linroot.find(".//*[@cameraNumber='"+str(self.camera)+"'][@ccdNumber='"+str(self.ccd)+"'][@ccdOutput='"+output+"']/linearityPoly")
		offsetx = float(poly.get('offsetx'))
		order = int(poly.get('order'))
		originx = float(poly.get('originx'))
		scalex = float(poly.get('scalex'))
		maxDomain = float(poly.get('maxDomain'))

		coeffs = []
		for elem in poly.findall('coeffs/coeff'):
			coeffs.append(float(elem.get("value")))

		p = np.polynomial.Polynomial(np.asarray(coeffs))

		img = img*p((img-originx)*scalex)*self.hdu[1].header['GAIN'+output]

		return img

# =========================================================================
# Correct for LDE undershoot
# =========================================================================

	def undershoot(self,img,output,verbose=True):
		if verbose:
			print("")
			print("Correcting for LDE undershoot")

		# Load undershoot model
		undertree = ET.parse(pkg_resources.resource_stream(__name__, 'data/tess-undershoot.xml'))
		underroot = undertree.getroot()

		undershoot = float(underroot.find(".//*[@cameraNumber='"+str(self.camera)+"'][@ccdNumber='"+str(self.ccd)+"'][@ccdOutput='"+output+"']").get('undershoot'))

		# Undershoot appears on the right
		if (output == "A") or (output == "C"):
			img[:,:,1:] = img[:,:,1:] + undershoot * img[:,:,:-1]
			img[:,:,0] = img[:,:,0]*(1+undershoot)

		# Undershoot appears on the left
		if (output == "B") or (output == "D"):
			img[:,:,:-1] = img[:,:,:-1] + undershoot * img[:,:,1:]
			img[:,:,-1] = img[:,:,-1]*(1+undershoot)

		return img

# =========================================================================
# Calculate smear
# =========================================================================

	def calcsmear(self,img,smrow,output,smear='alternate',smear_range=None,wings=None,verbose=True):
		if verbose:
			print("")
			if smear == 'alternate':
				print("Calculating photometric smear from TPF background")
			if smear == 'standard':
				print("Calculating photometric smear from collateral data")
			if smear == 'both':
				print("Calculating photometric smear using both TPF background and collateral data")

		#Origin and dimensions of the target pixel file relative to the smear data
		if output == 'A':
			sx1 = 44
		if output == 'B':
			sx1 = 556
		if output == 'C':
			sx1 = 1068
		if output == 'D':
			sx1 = 1580

		x1 = self.hdu[2].header['CRVAL1P']-1-sx1
		dx = self.hdu[2].header['NAXIS1']
		dy = self.hdu[2].header['NAXIS2']

		if smear_range == None or wings == None:
			print('Warning: Specify both smear_range and wings. Revert to calculating photometric smear from TPF background')
			smear = 'alternate'

		#Calculate the 'regular' smear correction
		smcor = np.nanmedian(smrow,axis=1)

		if smear == 'alternate' or smear == 'both':
			#Find which rows are below the bleed column
			medimg = np.nanmedian(img,axis=0)

			maxrows = []
			for col in np.arange(dx):
				maxrows.append(np.argmax(medimg[:,col][medimg[:,col] < 1e5]))

			maxrow = np.nanmin(maxrows) - 10 #10 row buffer to (hopefully) ensure we're always below the bleed column

			#Build a mask from the pixels with the lowest flux in each column
			smmask = np.zeros((dy,dx),dtype='bool')

			for col in np.arange(dx):
				idx = medimg[:maxrow,col] < np.percentile(medimg[:maxrow,col],20)
				smmask[:maxrow,col][idx] = True

			if smear == 'alternate':
				#From the background pixels calculate a smear correction + background for each column as a function of time
				smbkgd = np.nanmean(img.T[smmask.T].reshape(dx,int(np.sum(smmask)/dx),img.shape[0]),axis=1).T

				#Estimate the background and remove it from the smear correction. Here we can use the median level for the regular smear correction
				msk = (np.arange(smrow.shape[2]) < x1) | (np.arange(smrow.shape[2]) > x1 + dx)
				bkgd = np.nanmin(smbkgd,axis=1) - np.nanmedian(smcor[:,msk],axis=1)

				smcor[:,np.max([x1,0]):x1+dx] = smbkgd[:,-dx-np.min([x1,0]):] - np.tile(bkgd[:,np.newaxis],dx+np.min([x1,0]))
            
			if smear == 'both':
				left_smear, right_smear = smear_range[0] - (sx1 + x1) - 1, smear_range[1] - (sx1 + x1)
				left_wings, right_wings = wings[0] - (sx1 + x1) - 1, wings[1] - (sx1 + x1)

				#From the background pixels calculate a smear correction + background for each column as a function of time
				smbkgd = np.nanmean(img.T[smmask.T].reshape(dx,int(np.sum(smmask)/dx),img.shape[0]),axis=1).T

				#Estimate the background and remove it from the smear correction. Here we can use the median level for the regular smear correction
				msk = (np.arange(smrow.shape[2]) < x1) | (np.arange(smrow.shape[2]) > x1 + dx)
				bkgd = np.nanmin(smbkgd[:, left_wings: right_wings],axis=1) - np.nanmedian(smcor[:,msk],axis=1)

				smcor_copy = smbkgd[:,-dx-np.min([x1,0]):] - np.tile(bkgd[:,np.newaxis],dx+np.min([x1,0]))
				smcor[:,np.max([x1 + left_smear,0]):x1+right_smear] = smcor_copy[:,left_smear:right_smear]
		
        #Save calibration data for later use
		col1 = fits.Column(name=self.hdu[1].header['TTYPE1'], format=self.hdu[1].header['TFORM1'], unit=self.hdu[1].header['TUNIT1'], disp=self.hdu[1].header['TDISP1'], array=self.hdu[1].data['TIME'])
		col2 = fits.Column(name=self.hdu[1].header['TTYPE2'], format=self.hdu[1].header['TFORM2'], disp=self.hdu[1].header['TDISP2'], array=self.hdu[1].data['CADENCENO'])
		col3 = fits.Column(name='SMEAR', format=str(smcor.shape[1])+'D', unit='count', disp='F14.7', dim='('+str(smcor.shape[1])+')', array=smcor)

		cols = fits.ColDefs([col1,col2,col3])
		smearhdu = fits.BinTableHDU.from_columns(cols)

		path=self.path
		smearfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-smrow-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col-blf.fits'
		smearhdu.writeto(smearfile)		

		smrow[:,0,:] = smcor

		return smrow

# =========================================================================
# Remove smear
# =========================================================================
	
	def smear(self,img,sm,output,verbose=True):
		if verbose:
			print("")
			print("Removing photometric smear")

		#Origin and dimensions of the target pixel file relative to the smear data
		if output == 'A':
			sx1 = 44
		if output == 'B':
			sx1 = 556
		if output == 'C':
			sx1 = 1068
		if output == 'D':
			sx1 = 1580

		x1 = self.hdu[2].header['CRVAL1P']-1-sx1
		dx = self.hdu[2].header['NAXIS1']

		if sm.shape[0] == 2: #TPF goes over two outputs
			newsm = np.zeros((sm.shape[1],10,1024))
			newsm[:,:,:512] = sm[0]
			newsm[:,:,512:] = sm[1]
			sm = newsm
		else:
			sm = sm[0]		

		img = img - np.expand_dims(sm[:,x1:x1+dx],axis=1)

		return img

# =========================================================================
# Divide by flat field
# =========================================================================

	def flatfield(self,img,verbose=True):
		flatfile =  self.ddir+'tess2018323-'+str(self.camera)+'-'+str(self.ccd)+'-flat.fits'

		if verbose:
			print("")
			print("Loading flat field")
			print("... "+flatfile)
		
		flatfits = fits.open(flatfile)

		if verbose:
			print("Performing flatfield correction")

		#Define x,y and dimensions of flat field to cut out
		x1 = self.hdu[2].header['CRVAL1P']-1
		dx = self.hdu[2].header['NAXIS1']

		y1 = self.hdu[2].header['CRVAL2P']-1
		dy = self.hdu[2].header['NAXIS2']

		flat = flatfits[1].data[y1:y1+dy,x1:x1+dx]
		img = img/flat

		return img

# =========================================================================
# Convert to electrons per second
# =========================================================================

	def to_electrons_per_second(self,img,verbose=True):
		if verbose:
			print("")
			print("Converting to electrons per second")

		# Current units are electrons per exposure
		# Divide by seconds per exposure
		img = img/self.hdu[1].header['INT_TIME']

		return img
