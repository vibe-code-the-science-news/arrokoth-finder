#!/usr/bin/env python3
"""
Arrokoth (2014 MU69) Detection and Analysis Script
Analyzes HST FITS files to detect and track Arrokoth's motion

Usage:
    python arrokoth_finder.py --data_dir ./data --csv arrokoth.csv.csv
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
import warnings
warnings.filterwarnings('ignore')

class ArrokothFinder:
    def __init__(self, data_dir, csv_file):
        """
        Initialize the Arrokoth finder with data directory and observation catalog.
        
        Parameters:
        -----------
        data_dir : str
            Path to directory containing FITS files
        csv_file : str
            Path to MAST observation catalog CSV
        """
        self.data_dir = data_dir
        self.csv_file = csv_file
        self.observations = None
        self.arrokoth_obs = None
        self.load_catalog()
        
    def load_catalog(self):
        """Load and parse the MAST observation catalog."""
        print("Loading observation catalog...")
        
        # Read CSV, skipping comment lines
        with open(self.csv_file, 'r') as f:
            lines = f.readlines()
        
        # Find data start line
        data_start = 0
        for i, line in enumerate(lines):
            if not line.startswith('#') and 'dataproduct_type' in line:
                data_start = i
                break
        
        # Read actual data
        self.observations = pd.read_csv(self.csv_file, skiprows=data_start)
        
        # Filter for Arrokoth observations (1110113Y designation)
        self.arrokoth_obs = self.observations[
            self.observations['target_name'].str.contains('1110113Y', na=False)
        ].copy()
        
        # Convert times to datetime
        # MAST uses Modified Julian Date (MJD)
        self.arrokoth_obs['datetime'] = pd.to_datetime(
            (self.arrokoth_obs['t_min'] - 40587) * 86400, unit='s'
        )
        
        # Sort by observation time
        self.arrokoth_obs = self.arrokoth_obs.sort_values('t_min').reset_index(drop=True)
        
        print(f"Found {len(self.arrokoth_obs)} Arrokoth observations")
        print(f"Time range: {self.arrokoth_obs['datetime'].min()} to {self.arrokoth_obs['datetime'].max()}")
        
    def find_fits_file(self, obs_id):
        """
        Find FITS file for a given observation ID.
        
        Parameters:
        -----------
        obs_id : str
            HST observation ID (e.g., 'iciij8f2q')
            
        Returns:
        --------
        str or None : Path to FITS file if found
        """
        # Look for various FITS file types
        suffixes = ['_drz.fits', '_flt.fits', '_raw.fits', '.fits']
        
        for root, dirs, files in os.walk(self.data_dir):
            for suffix in suffixes:
                filename = f"{obs_id}{suffix}"
                if filename in files:
                    return os.path.join(root, filename)
        return None
    
    def load_image(self, fits_path):
        """
        Load HST FITS image and extract key information.
        
        Parameters:
        -----------
        fits_path : str
            Path to FITS file
            
        Returns:
        --------
        dict : Image data, WCS, and metadata
        """
        try:
            with fits.open(fits_path) as hdul:
                # Find the science extension (usually extension 1 for HST)
                sci_ext = 1 if len(hdul) > 1 else 0
                
                data = hdul[sci_ext].data
                header = hdul[sci_ext].header
                
                # Create WCS object for coordinate transformations
                wcs = WCS(header)
                
                return {
                    'data': data,
                    'header': header,
                    'wcs': wcs,
                    'exptime': header.get('EXPTIME', 370),
                    'filter': header.get('FILTER1', 'F350LP'),
                    'date_obs': header.get('DATE-OBS', 'Unknown')
                }
        except Exception as e:
            print(f"Error loading {fits_path}: {e}")
            return None
    
    def detect_sources(self, image_data, threshold_sigma=5):
        """
        Detect sources in the image using DAOStarFinder.
        
        Parameters:
        -----------
        image_data : numpy.ndarray
            2D image data
        threshold_sigma : float
            Detection threshold in sigma above background
            
        Returns:
        --------
        astropy.table.Table : Detected sources
        """
        # Calculate image statistics
        mean, median, std = sigma_clipped_stats(image_data, sigma=3.0)
        threshold = median + (threshold_sigma * std)
        
        # Detect sources
        daofind = DAOStarFinder(fwhm=3.0, threshold=threshold)
        sources = daofind(image_data - median)
        
        return sources
    
    def predict_arrokoth_position(self, wcs, obs_time, rough_ra=280.3, rough_dec=-20.9):
        """
        Predict Arrokoth's approximate pixel position in the image.
        This is a rough calculation - Arrokoth's orbit causes position changes.
        
        Parameters:
        -----------
        wcs : astropy.wcs.WCS
            World coordinate system of the image
        obs_time : float
            Observation time (MJD)
        rough_ra, rough_dec : float
            Approximate sky coordinates (degrees)
            
        Returns:
        --------
        tuple : (x, y) pixel coordinates
        """
        # This is a simplified prediction - in reality, you'd use orbital elements
        # For now, use the catalog coordinates as a starting point
        sky_coord = SkyCoord(ra=rough_ra*u.degree, dec=rough_dec*u.degree)
        x, y = wcs.world_to_pixel(sky_coord)
        return float(x), float(y)
    
    def analyze_observation(self, obs_row, show_plot=True):
        """
        Analyze a single observation to find and measure Arrokoth.
        
        Parameters:
        -----------
        obs_row : pandas.Series
            Single row from observations DataFrame
        show_plot : bool
            Whether to display diagnostic plots
            
        Returns:
        --------
        dict : Analysis results
        """
        obs_id = obs_row['obs_id']
        print(f"\nAnalyzing observation {obs_id}...")
        
        # Find FITS file
        fits_path = self.find_fits_file(obs_id)
        if not fits_path:
            print(f"FITS file not found for {obs_id}")
            return None
            
        print(f"Found FITS file: {fits_path}")
        
        # Load image
        img_data = self.load_image(fits_path)
        if not img_data:
            return None
            
        # Detect sources
        sources = self.detect_sources(img_data['data'])
        if sources is None:
            print("No sources detected")
            return None
            
        print(f"Detected {len(sources)} sources")
        
        # Predict Arrokoth position
        pred_x, pred_y = self.predict_arrokoth_position(
            img_data['wcs'], 
            obs_row['t_min'],
            obs_row['s_ra'], 
            obs_row['s_dec']
        )
        
        # Find sources near predicted position (within 50 pixels)
        distances = np.sqrt((sources['xcentroid'] - pred_x)**2 + 
                          (sources['ycentroid'] - pred_y)**2)
        nearby_sources = sources[distances < 50]
        
        if show_plot:
            self.plot_observation(img_data, sources, nearby_sources, pred_x, pred_y, obs_id)
        
        # Performance photometry on potential Arrokoth sources
        photometry_results = []
        for source in nearby_sources:
            aperture = CircularAperture((source['xcentroid'], source['ycentroid']), r=5.0)
            phot_table = aperture_photometry(img_data['data'], aperture)
            
            photometry_results.append({
                'x': source['xcentroid'],
                'y': source['ycentroid'], 
                'flux': phot_table['aperture_sum'][0],
                'mag': -2.5 * np.log10(phot_table['aperture_sum'][0] / img_data['exptime']),
                'distance_from_prediction': distances[sources == source][0]
            })
        
        return {
            'obs_id': obs_id,
            'fits_path': fits_path,
            'datetime': obs_row['datetime'],
            'predicted_pos': (pred_x, pred_y),
            'catalog_pos': (obs_row['s_ra'], obs_row['s_dec']),
            'total_sources': len(sources),
            'nearby_sources': len(nearby_sources),
            'photometry': photometry_results,
            'image_info': {
                'filter': img_data['filter'],
                'exptime': img_data['exptime'],
                'date_obs': img_data['date_obs']
            }
        }
    
    def plot_observation(self, img_data, all_sources, nearby_sources, pred_x, pred_y, obs_id):
        """Create diagnostic plot showing detected sources and predicted Arrokoth position."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Full image with all sources
        im1 = ax1.imshow(img_data['data'], origin='lower', cmap='gray', 
                        vmin=np.percentile(img_data['data'], 1),
                        vmax=np.percentile(img_data['data'], 99))
        ax1.scatter(all_sources['xcentroid'], all_sources['ycentroid'], 
                   c='red', s=20, alpha=0.7, marker='+')
        ax1.plot(pred_x, pred_y, 'bo', markersize=10, fillstyle='none', 
                label='Predicted Arrokoth')
        ax1.set_title(f'{obs_id} - Full Field ({len(all_sources)} sources)')
        ax1.legend()
        
        # Zoomed region around predicted position
        zoom_size = 100
        x_min = max(0, int(pred_x - zoom_size))
        x_max = min(img_data['data'].shape[1], int(pred_x + zoom_size))
        y_min = max(0, int(pred_y - zoom_size))
        y_max = min(img_data['data'].shape[0], int(pred_y + zoom_size))
        
        zoomed_data = img_data['data'][y_min:y_max, x_min:x_max]
        im2 = ax2.imshow(zoomed_data, origin='lower', cmap='gray',
                        vmin=np.percentile(zoomed_data, 1),
                        vmax=np.percentile(zoomed_data, 99))
        
        # Plot nearby sources in zoomed view
        for source in nearby_sources:
            x_zoom = source['xcentroid'] - x_min
            y_zoom = source['ycentroid'] - y_min
            ax2.plot(x_zoom, y_zoom, 'ro', markersize=8, fillstyle='none')
            ax2.text(x_zoom + 3, y_zoom + 3, f"{source['mag']:.1f}", 
                    color='red', fontsize=8)
        
        # Mark predicted position
        ax2.plot(pred_x - x_min, pred_y - y_min, 'bo', markersize=12, 
                fillstyle='none', linewidth=2, label='Predicted')
        ax2.set_title(f'Zoomed Region ({len(nearby_sources)} nearby sources)')
        ax2.legend()
        
        plt.tight_layout()
        plt.show()
    
    def analyze_all_observations(self, max_obs=10):
        """
        Analyze multiple Arrokoth observations to track its motion.
        
        Parameters:
        -----------
        max_obs : int
            Maximum number of observations to analyze
        """
        print(f"\n{'='*60}")
        print("ARROKOTH DETECTION AND TRACKING ANALYSIS")
        print(f"{'='*60}")
        
        results = []
        
        # Analyze a subset of observations
        obs_subset = self.arrokoth_obs.head(max_obs)
        
        for idx, obs in obs_subset.iterrows():
            result = self.analyze_observation(obs, show_plot=True)
            if result:
                results.append(result)
        
        # Summary analysis
        if results:
            self.plot_motion_summary(results)
            self.create_summary_report(results)
        
        return results
    
    def plot_motion_summary(self, results):
        """Plot Arrokoth's motion across multiple observations."""
        if len(results) < 2:
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Sky motion plot
        times = [r['datetime'] for r in results]
        ras = [r['catalog_pos'][0] for r in results]
        decs = [r['catalog_pos'][1] for r in results]
        
        ax1.plot(ras, decs, 'bo-', markersize=6)
        ax1.set_xlabel('RA (degrees)')
        ax1.set_ylabel('Dec (degrees)')
        ax1.set_title('Arrokoth Sky Motion')
        ax1.grid(True, alpha=0.3)
        
        # Add time labels
        for i, (ra, dec, time) in enumerate(zip(ras, decs, times)):
            ax1.annotate(f'{i+1}', (ra, dec), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8)
        
        # Detection statistics
        obs_ids = [r['obs_id'] for r in results]
        nearby_counts = [r['nearby_sources'] for r in results]
        
        ax2.bar(range(len(obs_ids)), nearby_counts)
        ax2.set_xlabel('Observation Number')
        ax2.set_ylabel('Nearby Sources Detected')
        ax2.set_title('Source Detection Statistics')
        ax2.set_xticks(range(len(obs_ids)))
        ax2.set_xticklabels([f"{i+1}" for i in range(len(obs_ids))], rotation=45)
        
        plt.tight_layout()
        plt.show()
    
    def create_summary_report(self, results):
        """Create a summary report of the analysis."""
        print(f"\n{'='*60}")
        print("ANALYSIS SUMMARY REPORT")
        print(f"{'='*60}")
        
        print(f"Observations analyzed: {len(results)}")
        print(f"Time span: {results[0]['datetime']} to {results[-1]['datetime']}")
        
        # Source detection statistics
        total_sources = sum(r['total_sources'] for r in results)
        nearby_sources = sum(r['nearby_sources'] for r in results)
        
        print(f"\nSource Detection:")
        print(f"  Total sources detected: {total_sources}")
        print(f"  Sources near predicted Arrokoth position: {nearby_sources}")
        print(f"  Average sources per image: {total_sources/len(results):.1f}")
        
        # Sky motion analysis
        if len(results) > 1:
            ra_range = max(r['catalog_pos'][0] for r in results) - min(r['catalog_pos'][0] for r in results)
            dec_range = max(r['catalog_pos'][1] for r in results) - min(r['catalog_pos'][1] for r in results)
            
            print(f"\nSky Motion:")
            print(f"  RA range: {ra_range*3600:.1f} arcseconds")
            print(f"  Dec range: {dec_range*3600:.1f} arcseconds")
            print(f"  Total motion: {np.sqrt(ra_range**2 + dec_range**2)*3600:.1f} arcseconds")
        
        print(f"\n{'='*60}")

def main():
    """Main function to run Arrokoth analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze HST observations to find Arrokoth')
    parser.add_argument('--data_dir', default='./data', help='Directory containing FITS files')
    parser.add_argument('--csv', default='arrokoth.csv.csv', help='MAST observation catalog CSV file')
    parser.add_argument('--max_obs', type=int, default=5, help='Maximum observations to analyze')
    parser.add_argument('--interactive', action='store_true', help='Interactive mode with plots')
    
    args = parser.parse_args()
    
    # Initialize finder
    finder = ArrokothFinder(args.data_dir, args.csv)
    
    # Run analysis
    results = finder.analyze_all_observations(max_obs=args.max_obs)
    
    print(f"\nAnalysis complete! Processed {len(results)} observations.")
    return results

if __name__ == "__main__":
    results = main()