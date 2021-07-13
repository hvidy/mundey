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

import warnings

class mundey_tpf(lightkurve.TessTargetPixelFile):


# =========================================================================
# Main routine to perform a new calibration of a TESS target pixel file
# =========================================================================

	def calibrate(self,verbose=True,smear='alternate',ddir='./',target=None,use_simbad=True):
		"""Performs a new calibration for a target pixel file
			Parameters
			----------
			smear : 'standard', or 'alternate'
				If the string 'standard' is passed, will calculate the smear correction directly from the collateral files.
				If the string 'alternate' is passed, smear correction will be estimated from the target pixel file.

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

		if verbose:
			print("Sector: "+str(self.sector))
			print("Camera: " +str(self.camera))
			print("CCD: " +str(self.ccd))
			print("Outputs: " +str(set(outputs)))
			print("")
			print("Loading collateral target pixel files")

		sm_cal = np.zeros((len(set(outputs)),flux.shape[0],10,512))
		tv_cal = np.zeros((len(set(outputs)),flux.shape[0],2078,11))

		for idx,output in enumerate(set(outputs)):

			smrowfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-smrow-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col.fits'
			if verbose:
				print("... loading smear row file "+smrowfile)

			smrow = fits.open(smrowfile)
			sm_cal[idx] = smrow[1].data['SMROW_RAW']

			tvcolfile = self.ddir+path.split('/')[-1].split('-')[0]+'-'+path.split('/')[-1].split('-')[1]+'-tvcol-'+str(self.camera)+'-'+str(self.ccd)+'-'+output.lower()+'-'+path.split('/')[-1].split('-')[3]+'-s_col.fits'
			if verbose:
				print("... loading trailing virtual column file "+tvcolfile)

			tvcol = fits.open(tvcolfile)
			tv_cal[idx] = tvcol[1].data['TVCOL_RAW']


		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux[:,:,indices] = self.fixedoffset_meanblack(flux[:,:,indices],output,verbose=verbose)
			sm_cal[idx] = self.fixedoffset_meanblack(sm_cal[idx],output,verbose=False)
			tv_cal[idx] = self.fixedoffset_meanblack(tv_cal[idx],output,verbose=False)
	

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
			sm_cal[idx] = self.twodblack(sm_cal[idx],black2dimg[2058:2068,sx1:sx2],verbose=False)
			tv_cal[idx] = self.twodblack(tv_cal[idx],black2dimg[:,tx1:tx2],verbose=False)
		
		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux,sm_cal[idx] = self.onedblack(flux,sm_cal[idx],tv_cal[idx],verbose=verbose)

		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]
		
			flux[:,:,indices] = self.linearity_gain(flux[:,:,indices],output,verbose=verbose)
			sm_cal[idx] = self.linearity_gain(sm_cal[idx],output,verbose=False)

		for idx,output in enumerate(set(outputs)):
			indices = [i for i, x in enumerate(outputs) if x == output]

			flux[:,:,indices] = self.undershoot(flux[:,:,indices],output,verbose=verbose)
			sm_cal[idx] = self.undershoot(sm_cal[idx],output,verbose=False)


		flux = self.smear(flux,sm_cal,outputs[0],smear=smear,verbose=verbose)

		flux = self.flatfield(flux,verbose=verbose)

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
# Calculate 1D black and remove it
# =========================================================================

	def onedblack(self,img,smrow,tvcol,verbose=True):
		if verbose:
			print("")
			print("Performing 1-D black correction")
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

		for frame in tqdm(np.arange(tvcol.shape[0]),desc='Frame'):
			aics = []
			y = np.nanmean(tvcol[frame,:,:-1],axis=1) #Average over all but the last column, which seems to be off

			for k in ks:
				p = np.polynomial.Polynomial.fit(x,y,k)
				aic = 2*k + nrows*np.log(np.sum((y-p(x))**2))
				aics.append(aic)

			order = ks[np.argmin(np.asarray(aics))]
			p = np.polynomial.Polynomial.fit(x,y,order)

			black1d = p(x)

			img[frame,:,:] = img[frame,:,:] - np.expand_dims(black1d[y1:y1+dy],axis=1)
			smrow[frame,:,:] = smrow[frame,:,:] - np.expand_dims(black1d[2058:2068],axis=1)

		return img,smrow

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
# Calculate and correct smear
# =========================================================================

	def smear(self,img,smrow,output,smear='alternate',verbose=True):
		if verbose:
			print("")
			print("Correcting for photometric smear")

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

		if smrow.shape[0] == 2: #TPF goes over two outputs
			newsmrow = np.zeros((smrow.shape[1],10,1024))
			newsmrow[:,:,:512] = smrow[0]
			newsmrow[:,:,512] = smrow[1]
			smrow = newsmrow
		else:
			smrow = smrow[0]

		#Calculate the 'regular' smear correction
		smcor = np.nanmedian(smrow,axis=1)

		if smear == 'alternate':
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

			#From the background pixels calculate a smear correction + background for each column as a function of time
			smbkgd = np.nanmean(img.T[smmask.T].reshape(dx,int(np.sum(smmask)/dx),img.shape[0]),axis=1).T

			#Estimate the background and remove it from the smear correction. Here we can use the median level for the regular smear correction
			msk = (np.arange(smrow.shape[2]) < x1) | (np.arange(smrow.shape[2]) > x1 + dx)
			bkgd = np.nanmin(smbkgd,axis=1) - np.nanmedian(smcor[:,msk],axis=1)
			sm = smbkgd - np.tile(bkgd[:,np.newaxis],dx)


		if smear == 'standard':
			sm = smcor[:,x1:x1+dx]

		#Apply the new correction
		img = img - np.expand_dims(sm,axis=1)


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