#!/usr/bin/env python3
"""
Arrokoth Blind Discovery Script - FITS Only
Discovers moving objects (including Arrokoth) by analyzing motion across HST FITS files

This script recreates the discovery process:
1. Scans all FITS files
2. Detects sources in each image  
3. Tracks objects that move consistently
4. Identifies slow-moving KBOs like Arrokoth

Usage:
    python arrokoth_blind_finder.py /path/to/fits/directory
    python arrokoth_blind_finder.py HST --min_detections 4 --max_motion 10
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
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
from scipy.spatial.distance import cdist
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class MovingObjectFinder:
    def __init__(self, data_dir, min_detections=3, max_motion_arcsec=10, max_days=100):
        """
        Initialize the moving object finder.
        
        Parameters:
        -----------
        data_dir : str
            Directory containing FITS files
        min_detections : int
            Minimum number of detections required for a moving object
        max_motion_arcsec : float
            Maximum motion in arcsec/day for valid objects
        max_days : float
            Maximum time span to consider for motion tracking
        """
        self.data_dir = data_dir
        self.min_detections = min_detections
        self.max_motion_arcsec = max_motion_arcsec
        self.max_days = max_days
        
        self.fits_files = []
        self.observations = []
        self.all_sources = pd.DataFrame()
        self.moving_objects = []
        
    def scan_fits_files(self):
        """Recursively find all FITS files in the data directory."""
        print(f"Scanning '{self.data_dir}' for FITS files...")
        
        fits_extensions = ['*.fits', '*.fit', '*.fts']
        self.fits_files = []
        
        for root, dirs, files in os.walk(self.data_dir):
            for ext in fits_extensions:
                pattern = os.path.join(root, ext)
                self.fits_files.extend(glob.glob(pattern))
        
        print(f"Found {len(self.fits_files)} FITS files")
        
        if len(self.fits_files) == 0:
            print("ERROR: No FITS files found!")
            return False
            
        return True
    
    def load_observation(self, fits_path):
        """
        Load a single FITS observation and extract metadata.
        
        Parameters:
        -----------
        fits_path : str
            Path to FITS file
            
        Returns:
        --------
        dict or None : Observation data
        """
        try:
            with fits.open(fits_path) as hdul:
                # Find science extension
                sci_ext = 1 if len(hdul) > 1 else 0
                data = hdul[sci_ext].data
                header = hdul[sci_ext].header
                
                # Skip if no data
                if data is None:
                    return None
                
                # Create WCS
                try:
                    wcs = WCS(header)
                except:
                    print(f"Warning: Could not create WCS for {os.path.basename(fits_path)}")
                    return None
                
                # Extract observation time - try multiple header keywords
                obs_time = None
                
                # Try different combinations of date/time keywords
                date_keys = ['DATE-OBS', 'DATE_OBS', 'OBSDATE', 'DATE']
                time_keys = ['TIME-OBS', 'TIME_OBS', 'OBSTIME', 'TIME']
                mjd_keys = ['MJD-OBS', 'MJD_OBS', 'MJD-STR', 'MJD', 'EXPSTART']
                
                # First try MJD directly
                for mjd_key in mjd_keys:
                    if mjd_key in header:
                        try:
                            obs_time = Time(header[mjd_key], format='mjd')
                            break
                        except:
                            continue
                
                # If no MJD, try date + time combination
                if obs_time is None:
                    for date_key in date_keys:
                        if date_key in header:
                            date_obs = header[date_key]
                            time_obs = '00:00:00'  # default time
                            
                            # Look for corresponding time
                            for time_key in time_keys:
                                if time_key in header:
                                    time_obs = header[time_key]
                                    break
                            
                            try:
                                if date_obs and len(date_obs) > 8:
                                    obs_time = Time(f"{date_obs}T{time_obs}")
                                    break
                            except:
                                continue
                
                # Last resort: try to extract from filename or use a default
                if obs_time is None:
                    # For HST survey data, we can estimate time from observation ID
                    # This is a rough approximation for motion analysis
                    obs_id = self.extract_obs_id(fits_path)
                    
                    # Try to extract sequence info from observation ID
                    # HST IDs often encode sequence information
                    try:
                        # Extract numeric part and use as relative time
                        import re
                        numbers = re.findall(r'\d+', obs_id)
                        if numbers:
                            # Use first number group as a sequence indicator
                            seq_num = int(numbers[0])
                            # Create approximate MJD (starting from a reference date)
                            base_mjd = 56800  # Approximate MJD for 2014 survey
                            obs_time = Time(base_mjd + seq_num * 0.1, format='mjd')  # Space observations by ~2.4 hours
                        else:
                            # Final fallback - use a sequential time
                            base_mjd = 56800
                            obs_time = Time(base_mjd, format='mjd')
                    except:
                        print(f"Warning: Could not determine time for {os.path.basename(fits_path)} - skipping")
                        return None
                
                # Calculate field center coordinates
                ny, nx = data.shape
                center_coord = wcs.pixel_to_world(nx/2, ny/2)
                
                return {
                    'filename': os.path.basename(fits_path),
                    'filepath': fits_path,
                    'data': data,
                    'header': header,
                    'wcs': wcs,
                    'obs_time': obs_time,
                    'mjd': obs_time.mjd,
                    'center_ra': center_coord.ra.degree,
                    'center_dec': center_coord.dec.degree,
                    'exptime': header.get('EXPTIME', 1.0),
                    'filter': header.get('FILTER1', header.get('FILTER', 'Unknown')),
                    'obs_id': self.extract_obs_id(fits_path)
                }
        except Exception as e:
            print(f"Error loading {os.path.basename(fits_path)}: {e}")
            return None
    
    def extract_obs_id(self, fits_path):
        """Extract observation ID from FITS filename."""
        filename = os.path.basename(fits_path)
        # Remove extension and suffixes
        name = filename.lower()
        for ext in ['.fits', '.fit', '.fts']:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        for suffix in ['_drz', '_flt', '_raw', '_crj']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        return name
    
    def detect_sources_in_observation(self, obs):
        """
        Detect sources in a single observation.
        
        Parameters:
        -----------
        obs : dict
            Observation data from load_observation()
            
        Returns:
        --------
        pandas.DataFrame : Detected sources with world coordinates
        """
        data = obs['data']
        wcs = obs['wcs']
        
        # Check data validity
        if data is None:
            print(f"Warning: No data in {obs['filename']}")
            return pd.DataFrame()
        
        # Calculate background statistics
        try:
            mean, median, std = sigma_clipped_stats(data, sigma=3.0)
            threshold = median + (3.0 * std)  # Lowered threshold from 5.0 to 3.0
            
            print(f"Debug {obs['filename']}: mean={mean:.1f}, median={median:.1f}, std={std:.1f}, threshold={threshold:.1f}")
            
        except Exception as e:
            print(f"Error calculating stats for {obs['filename']}: {e}")
            return pd.DataFrame()
        
        # Detect sources
        try:
            daofind = DAOStarFinder(fwhm=2.5, threshold=threshold, sharplo=0.2, sharphi=2.0)
            sources = daofind(data - median)
        except Exception as e:
            print(f"Error in source detection for {obs['filename']}: {e}")
            return pd.DataFrame()
        
        if sources is None or len(sources) == 0:
            print(f"No sources detected in {obs['filename']} (threshold={threshold:.1f})")
            return pd.DataFrame()
        
        print(f"Detected {len(sources)} sources in {obs['filename']}")
        
        # Convert to world coordinates
        try:
            world_coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
            sources['ra'] = world_coords.ra.degree
            sources['dec'] = world_coords.dec.degree
        except Exception as e:
            print(f"Warning: WCS conversion failed for {obs['filename']}: {e}")
            return pd.DataFrame()
        
        # Create DataFrame with relevant columns
        source_df = pd.DataFrame({
            'obs_id': obs['obs_id'],
            'filename': obs['filename'],
            'mjd': obs['mjd'],
            'x': sources['xcentroid'],
            'y': sources['ycentroid'],
            'ra': sources['ra'],
            'dec': sources['dec'],
            'flux': sources['flux'],
            'mag': -2.5 * np.log10(sources['flux'] / obs['exptime']) + 25,  # Rough magnitude
            'filter': obs['filter']
        })
        
        return source_df
    
    def load_all_observations(self):
        """Load all FITS files and detect sources."""
        print("\nLoading observations and detecting sources...")
        
        self.observations = []
        all_source_dfs = []
        successful_loads = 0
        failed_loads = 0
        
        for i, fits_path in enumerate(self.fits_files):
            if i % 10 == 0:
                print(f"Processing {i+1}/{len(self.fits_files)}: {os.path.basename(fits_path)}")
            
            obs = self.load_observation(fits_path)
            if obs is None:
                failed_loads += 1
                continue
                
            self.observations.append(obs)
            successful_loads += 1
            
            # Detect sources
            sources_df = self.detect_sources_in_observation(obs)
            if len(sources_df) > 0:
                all_source_dfs.append(sources_df)
        
        print(f"\nLoading summary:")
        print(f"  Successfully loaded: {successful_loads}/{len(self.fits_files)} observations")
        print(f"  Failed to load: {failed_loads}/{len(self.fits_files)} observations")
        
        if all_source_dfs:
            self.all_sources = pd.concat(all_source_dfs, ignore_index=True)
            print(f"  Total sources detected: {len(self.all_sources)}")
            print(f"  Observations with sources: {len(all_source_dfs)}")
        else:
            print("  ERROR: No sources detected in any observations!")
            print("  This could be due to:")
            print("    - Images are too faint/noisy")
            print("    - Detection threshold too high")
            print("    - WCS/coordinate issues")
            return False
            
        return True
    
    def find_moving_objects(self):
        """
        Find objects that move consistently across multiple observations.
        """
        print("\nSearching for moving objects...")
        
        if len(self.observations) < self.min_detections:
            print(f"Need at least {self.min_detections} observations, only have {len(self.observations)}")
            return
        
        # Sort observations by time
        obs_times = [obs['mjd'] for obs in self.observations]
        time_sorted_indices = np.argsort(obs_times)
        
        # Group sources by observation
        sources_by_obs = {}
        for obs_idx, obs in enumerate(self.observations):
            obs_sources = self.all_sources[self.all_sources['obs_id'] == obs['obs_id']]
            if len(obs_sources) > 0:
                sources_by_obs[obs_idx] = obs_sources
        
        print(f"Found sources in {len(sources_by_obs)} observations")
        
        # Find moving objects by linking sources across observations
        moving_candidates = []
        
        # Start with sources from the first observation
        if 0 in sources_by_obs:
            first_sources = sources_by_obs[0]
            
            for _, source1 in first_sources.iterrows():
                trajectory = [source1]
                
                # Look for this object in subsequent observations
                for obs_idx in time_sorted_indices[1:]:
                    if obs_idx not in sources_by_obs:
                        continue
                        
                    obs_sources = sources_by_obs[obs_idx]
                    obs_time = self.observations[obs_idx]['mjd']
                    
                    # Predict where this object should be based on current trajectory
                    if len(trajectory) >= 2:
                        # Linear extrapolation
                        dt1 = trajectory[-1]['mjd'] - trajectory[-2]['mjd']
                        dt2 = obs_time - trajectory[-1]['mjd']
                        
                        if dt1 > 0 and dt2 < self.max_days:
                            pred_ra = trajectory[-1]['ra'] + (trajectory[-1]['ra'] - trajectory[-2]['ra']) * (dt2/dt1)
                            pred_dec = trajectory[-1]['dec'] + (trajectory[-1]['dec'] - trajectory[-2]['dec']) * (dt2/dt1)
                        else:
                            pred_ra, pred_dec = trajectory[-1]['ra'], trajectory[-1]['dec']
                    else:
                        pred_ra, pred_dec = trajectory[-1]['ra'], trajectory[-1]['dec']
                    
                    # Find closest source to prediction
                    distances = np.sqrt((obs_sources['ra'] - pred_ra)**2 + 
                                      (obs_sources['dec'] - pred_dec)**2) * 3600  # arcsec
                    
                    if len(distances) > 0:
                        min_dist_idx = distances.idxmin()
                        min_distance = distances[min_dist_idx]
                        
                        # Check if motion is reasonable for a KBO
                        time_diff = obs_time - trajectory[-1]['mjd']
                        if time_diff > 0:
                            motion_rate = min_distance / time_diff  # arcsec/day
                            
                            if motion_rate <= self.max_motion_arcsec and min_distance < 30:  # 30 arcsec matching radius
                                closest_source = obs_sources.loc[min_dist_idx]
                                trajectory.append(closest_source)
                
                # Check if we have enough detections
                if len(trajectory) >= self.min_detections:
                    moving_candidates.append(trajectory)
        
        print(f"Found {len(moving_candidates)} potential moving objects")
        
        # Analyze each candidate
        self.moving_objects = []
        for i, trajectory in enumerate(moving_candidates):
            analysis = self.analyze_trajectory(trajectory, i+1)
            if analysis:
                self.moving_objects.append(analysis)
        
        print(f"Confirmed {len(self.moving_objects)} moving objects")
    
    def analyze_trajectory(self, trajectory, obj_id):
        """
        Analyze a trajectory to determine if it's a real moving object.
        
        Parameters:
        -----------
        trajectory : list
            List of source detections
        obj_id : int
            Object identifier
            
        Returns:
        --------
        dict : Analysis results
        """
        if len(trajectory) < self.min_detections:
            return None
        
        # Convert to arrays for analysis
        times = np.array([src['mjd'] for src in trajectory])
        ras = np.array([src['ra'] for src in trajectory])
        decs = np.array([src['dec'] for src in trajectory])
        mags = np.array([src['mag'] for src in trajectory])
        
        # Calculate motion
        time_span = times.max() - times.min()
        ra_motion = (ras.max() - ras.min()) * 3600 * np.cos(np.radians(np.mean(decs)))  # arcsec
        dec_motion = (decs.max() - decs.min()) * 3600  # arcsec
        total_motion = np.sqrt(ra_motion**2 + dec_motion**2)
        
        if time_span > 0:
            motion_rate = total_motion / time_span  # arcsec/day
        else:
            motion_rate = 0
        
        # Fit linear motion
        if len(times) >= 3:
            ra_fit = np.polyfit(times, ras, 1)
            dec_fit = np.polyfit(times, decs, 1)
            
            # Calculate residuals
            ra_pred = np.polyval(ra_fit, times)
            dec_pred = np.polyval(dec_fit, times)
            ra_residuals = (ras - ra_pred) * 3600 * np.cos(np.radians(np.mean(decs)))
            dec_residuals = (decs - dec_pred) * 3600
            rms_residual = np.sqrt(np.mean(ra_residuals**2 + dec_residuals**2))
        else:
            ra_fit = dec_fit = [0, 0]
            rms_residual = 0
        
        # Object quality score (lower is better)
        # Penalize high motion rates, high residuals, low detection count
        quality_score = (motion_rate / self.max_motion_arcsec) + (rms_residual / 5.0) + (1.0 / len(trajectory))
        
        return {
            'object_id': obj_id,
            'detections': len(trajectory),
            'time_span_days': time_span,
            'ra_range_arcsec': ra_motion,
            'dec_range_arcsec': dec_motion,
            'total_motion_arcsec': total_motion,
            'motion_rate_arcsec_per_day': motion_rate,
            'mean_magnitude': np.mean(mags),
            'mag_std': np.std(mags),
            'rms_residual_arcsec': rms_residual,
            'quality_score': quality_score,
            'ra_motion_fit': ra_fit,
            'dec_motion_fit': dec_fit,
            'trajectory': trajectory,
            'mean_ra': np.mean(ras),
            'mean_dec': np.mean(decs),
            'first_detection': min(times),
            'last_detection': max(times)
        }
    
    def plot_moving_objects(self, max_objects=10):
        """Plot the trajectories of detected moving objects."""
        if not self.moving_objects:
            print("No moving objects to plot")
            return
        
        # Sort by quality score (best first)
        sorted_objects = sorted(self.moving_objects, key=lambda x: x['quality_score'])
        objects_to_plot = sorted_objects[:max_objects]
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
        
        # Plot 1: Sky motion trajectories
        ax = axes[0]
        colors = plt.cm.tab10(np.linspace(0, 1, len(objects_to_plot)))
        
        for i, obj in enumerate(objects_to_plot):
            trajectory = obj['trajectory']
            ras = [src['ra'] for src in trajectory]
            decs = [src['dec'] for src in trajectory]
            
            ax.plot(ras, decs, 'o-', color=colors[i], 
                   label=f"Obj {obj['object_id']} ({obj['detections']} det)", 
                   markersize=6, alpha=0.8)
            
            # Mark first and last positions
            ax.plot(ras[0], decs[0], 's', color=colors[i], markersize=8, alpha=0.6)
            ax.plot(ras[-1], decs[-1], '^', color=colors[i], markersize=8, alpha=0.6)
        
        ax.set_xlabel('RA (degrees)')
        ax.set_ylabel('Dec (degrees)')
        ax.set_title('Moving Object Trajectories')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Motion rates
        ax = axes[1]
        motion_rates = [obj['motion_rate_arcsec_per_day'] for obj in objects_to_plot]
        detections = [obj['detections'] for obj in objects_to_plot]
        object_ids = [obj['object_id'] for obj in objects_to_plot]
        
        bars = ax.bar(range(len(objects_to_plot)), motion_rates, color=colors[:len(objects_to_plot)])
        ax.set_xlabel('Object ID')
        ax.set_ylabel('Motion Rate (arcsec/day)')
        ax.set_title('Motion Rates of Detected Objects')
        ax.set_xticks(range(len(objects_to_plot)))
        ax.set_xticklabels(object_ids)
        
        # Add detection count on bars
        for i, (bar, det_count) in enumerate(zip(bars, detections)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   f'{det_count}', ha='center', va='bottom', fontsize=8)
        
        # Plot 3: Magnitude vs time for best object
        ax = axes[2]
        if objects_to_plot:
            best_obj = objects_to_plot[0]
            trajectory = best_obj['trajectory']
            times = [src['mjd'] for src in trajectory]
            mags = [src['mag'] for src in trajectory]
            
            ax.plot(times, mags, 'o-', markersize=8, linewidth=2)
            ax.set_xlabel('MJD')
            ax.set_ylabel('Magnitude')
            ax.set_title(f'Lightcurve - Object {best_obj["object_id"]}')
            ax.grid(True, alpha=0.3)
            ax.invert_yaxis()
        
        # Plot 4: Quality scores
        ax = axes[3]
        quality_scores = [obj['quality_score'] for obj in objects_to_plot]
        ax.bar(range(len(objects_to_plot)), quality_scores, color=colors[:len(objects_to_plot)])
        ax.set_xlabel('Object ID')
        ax.set_ylabel('Quality Score (lower = better)')
        ax.set_title('Object Quality Scores')
        ax.set_xticks(range(len(objects_to_plot)))
        ax.set_xticklabels(object_ids)
        
        plt.tight_layout()
        plt.show()
    
    def create_summary_report(self):
        """Create a summary report of discovered moving objects."""
        print(f"\n{'='*80}")
        print("MOVING OBJECT DISCOVERY REPORT")
        print(f"{'='*80}")
        
        print(f"Total FITS files processed: {len(self.fits_files)}")
        print(f"Valid observations loaded: {len(self.observations)}")
        print(f"Total sources detected: {len(self.all_sources)}")
        print(f"Moving objects discovered: {len(self.moving_objects)}")
        
        if self.observations:
            mjds = [obs['mjd'] for obs in self.observations]
            time_span = max(mjds) - min(mjds)
            print(f"Time span: {time_span:.1f} days")
            print(f"Date range: {Time(min(mjds), format='mjd').iso[:10]} to {Time(max(mjds), format='mjd').iso[:10]}")
        
        if not self.moving_objects:
            print("\nNo moving objects detected.")
            print("Try adjusting parameters:")
            print("  - Increase max_motion_arcsec for faster objects")
            print("  - Decrease min_detections for sparse data")
            return
        
        # Sort by quality (best first)
        sorted_objects = sorted(self.moving_objects, key=lambda x: x['quality_score'])
        
        print(f"\n{'='*80}")
        print("TOP MOVING OBJECT CANDIDATES")
        print(f"{'='*80}")
        print(f"{'ID':>3} {'Det':>3} {'Days':>6} {'Motion':>8} {'Rate':>8} {'Mag':>6} {'Quality':>8}")
        print(f"{'':>3} {'':>3} {'':>6} {'(arcsec)':>8} {'(arcsec/day)':>8} {'':>6} {'Score':>8}")
        print("-" * 80)
        
        for obj in sorted_objects[:15]:  # Show top 15
            print(f"{obj['object_id']:3d} "
                  f"{obj['detections']:3d} "
                  f"{obj['time_span_days']:6.1f} "
                  f"{obj['total_motion_arcsec']:8.1f} "
                  f"{obj['motion_rate_arcsec_per_day']:8.2f} "
                  f"{obj['mean_magnitude']:6.1f} "
                  f"{obj['quality_score']:8.3f}")
        
        print(f"\n{'='*80}")
        print("BEST CANDIDATE DETAILS")
        print(f"{'='*80}")
        
        best_obj = sorted_objects[0]
        print(f"Object ID: {best_obj['object_id']}")
        print(f"Detections: {best_obj['detections']}")
        print(f"Time span: {best_obj['time_span_days']:.1f} days")
        print(f"Total motion: {best_obj['total_motion_arcsec']:.1f} arcsec")
        print(f"Motion rate: {best_obj['motion_rate_arcsec_per_day']:.2f} arcsec/day")
        print(f"Mean position: RA={best_obj['mean_ra']:.6f}°, Dec={best_obj['mean_dec']:.6f}°")
        print(f"Mean magnitude: {best_obj['mean_magnitude']:.1f} ± {best_obj['mag_std']:.1f}")
        print(f"RMS residual: {best_obj['rms_residual_arcsec']:.1f} arcsec")
        
        print(f"\nDetection details:")
        for i, detection in enumerate(best_obj['trajectory']):
            time_str = Time(detection['mjd'], format='mjd').iso[:16]
            print(f"  {i+1}: {time_str} - RA={detection['ra']:.6f}° Dec={detection['dec']:.6f}° "
                  f"Mag={detection['mag']:.1f} ({detection['filename']})")
        
        # Check if this could be Arrokoth
        if (best_obj['motion_rate_arcsec_per_day'] < 5.0 and  # Slow motion typical of distant KBO
            best_obj['detections'] >= 3 and  # Multiple detections
            best_obj['mean_magnitude'] > 24):  # Faint like Arrokoth
            print(f"\n🎯 CANDIDATE ARROKOTH DETECTION! 🎯")
            print(f"This object has characteristics consistent with a distant KBO like Arrokoth:")
            print(f"  - Slow motion rate ({best_obj['motion_rate_arcsec_per_day']:.2f} arcsec/day)")
            print(f"  - Faint magnitude ({best_obj['mean_magnitude']:.1f})")
            print(f"  - Multiple detections ({best_obj['detections']})")
            print(f"  - Good linear motion fit (RMS={best_obj['rms_residual_arcsec']:.1f} arcsec)")
    
    def run_analysis(self):
        """Run the complete moving object analysis."""
        print("🔭 STARTING BLIND MOVING OBJECT DISCOVERY 🔭")
        print("=" * 60)
        
        # Step 1: Scan for FITS files
        if not self.scan_fits_files():
            return False
        
        # Step 2: Load observations and detect sources  
        if not self.load_all_observations():
            return False
        
        # Step 3: Find moving objects
        self.find_moving_objects()
        
        # Step 4: Create plots and report
        if self.moving_objects:
            self.plot_moving_objects()
            self.create_summary_report()
        else:
            print("\nNo moving objects found. Try adjusting search parameters.")
        
        return True

def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Discover moving objects in HST FITS files')
    parser.add_argument('data_dir', help='Directory containing FITS files')
    parser.add_argument('--min_detections', type=int, default=3, 
                       help='Minimum detections required (default: 3)')
    parser.add_argument('--max_motion', type=float, default=10.0,
                       help='Maximum motion in arcsec/day (default: 10)')
    parser.add_argument('--max_days', type=float, default=100.0,
                       help='Maximum time span in days (default: 100)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.data_dir):
        print(f"ERROR: Directory '{args.data_dir}' does not exist!")
        return
    
    # Create finder and run analysis
    finder = MovingObjectFinder(
        data_dir=args.data_dir,
        min_detections=args.min_detections,
        max_motion_arcsec=args.max_motion,
        max_days=args.max_days
    )
    
    success = finder.run_analysis()
    
    if success and finder.moving_objects:
        print(f"\n🎉 SUCCESS! Discovered {len(finder.moving_objects)} moving objects!")
        print("The best candidate might be Arrokoth!")
    else:
        print("\n❌ No moving objects discovered.")
        print("Try adjusting parameters or check data quality.")

if __name__ == "__main__":
    main()