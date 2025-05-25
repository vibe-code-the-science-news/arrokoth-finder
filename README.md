# 🚀 Arrokoth Discovery Project 🔭

_Recreating the discovery of the most distant object ever explored by a spacecraft_

![Arrokoth Header Image](https://newhorizons.jhuapl.edu/Galleries/Featured-Images/pics800wide/CA06_color_m-h_desmear_destrip_contrast-selected%20cover.png)
_Arrokoth (2014 MU69) as seen by New Horizons on January 1, 2019 - NASA/Johns Hopkins APL/Southwest Research Institute_

## 🌟 What is Arrokoth?

**Arrokoth** (formerly nicknamed "Ultima Thule") is a small, ancient object in the Kuiper Belt that made history on January 1, 2019, when NASA's New Horizons spacecraft conducted the most distant flyby in the history of space exploration. Located over 4 billion miles from Earth, this contact binary object revealed secrets about the early formation of our solar system.

But before New Horizons could visit Arrokoth, **it had to be discovered first**.

## 🔍 The Original Discovery Story

In 2014, with New Horizons speeding toward Pluto, mission planners faced a challenge: they wanted to visit another object after Pluto, but needed to find one that the spacecraft could actually reach with its limited remaining fuel.

**The Search Campaign:**

- **194 Hubble Space Telescope orbits** systematically scanned the sky
- **Marc Buie and the New Horizons team** analyzed thousands of images
- **5 potential targets** were discovered in the vast darkness
- **One faint moving dot** became 2014 MU69 - later named Arrokoth

This was **detective work at its finest** - finding a 22-mile-wide object 4 billion miles away by watching it slowly drift against the background stars.

## 🎯 What This Project Does

This project **recreates the original discovery process** using the actual Hubble Space Telescope data that found Arrokoth. We've built two complementary approaches:

### 1. **Guided Discovery Script** (`arrokoth_finder.py`)

- Uses the MAST observation catalog to focus on known Arrokoth observations
- Validates the discovery by confirming Arrokoth's position and motion
- Perfect for understanding how targeted follow-up observations work

### 2. **Blind Discovery Script** (`arrokoth_blind_finder.py`) ⭐

- **Recreates the original discovery from scratch**
- Analyzes HST survey images without knowing where Arrokoth should be
- Automatically detects **all moving objects** in the field
- **Rediscovers Arrokoth** based purely on its motion signature
- This is the **real deal** - digital archaeology of astronomical discovery!

## 🛠️ How It Works

### The Science Behind Moving Object Detection

1. **Source Detection**: Find all point sources (stars, galaxies, asteroids) in each image
2. **Coordinate Registration**: Convert pixel positions to sky coordinates using astrometry
3. **Motion Tracking**: Link sources across multiple observations taken over days/weeks
4. **Orbit Analysis**: Identify objects with consistent, slow motion typical of distant Kuiper Belt Objects
5. **Candidate Ranking**: Score objects based on detection quality, motion consistency, and brightness

### Key Characteristics That Identify Arrokoth:

- **Slow Motion**: ~2-5 arcseconds per day (much slower than nearby asteroids)
- **Faint Brightness**: Magnitude ~26-27 (near the detection limit)
- **Linear Trajectory**: Consistent motion in a straight line
- **Multiple Detections**: Appears reliably across many observations

## 📊 What You'll Discover

Running the blind discovery script on real HST data, you should find:

- **Several moving objects** in the survey field
- **Background stars** that don't move (reference frame)
- **Arrokoth itself** as the highest-quality slow-moving candidate
- **Motion plots** showing each object's path across the sky
- **Detailed analysis** of discovery statistics and object properties

## 🚀 Quick Start

### Prerequisites

```bash
pip install astropy photutils pandas matplotlib numpy scipy
```

### Basic Usage

```bash
# Recreate the blind discovery (most exciting!)
python arrokoth_blind_finder.py /path/to/HST/data

# Validate known Arrokoth observations
python arrokoth_finder.py /path/to/HST/data --csv arrokoth_observations.csv
```

### Sample Output

```
🔭 STARTING BLIND MOVING OBJECT DISCOVERY 🔭
Found 115 FITS files
Detected 12,847 total sources
Moving objects discovered: 4

🎯 CANDIDATE ARROKOTH DETECTION! 🎯
Object ID: 1
Detections: 8
Motion rate: 3.2 arcsec/day
Mean magnitude: 26.1
```

## 📈 Scientific Impact

This project demonstrates:

- **Modern astronomical survey techniques** for discovering small solar system bodies
- **The power of time-domain astronomy** - finding objects by their motion
- **Data mining approaches** for extracting discoveries from large datasets
- **The incredible precision** required for spacecraft navigation to distant targets

## 🎓 Educational Value

Perfect for:

- **Astronomy students** learning about observational techniques
- **Data science applications** in astronomical research
- **Understanding spacecraft mission planning** and target selection
- **Appreciating the technical challenges** of deep space exploration

## 🏆 The Bigger Picture

Arrokoth's discovery represents a pinnacle of **human curiosity and technical achievement**:

1. **Building Hubble Space Telescope** - decades of engineering
2. **Launching New Horizons** - precision trajectory planning
3. **Discovering the target** - finding a needle in a cosmic haystack
4. **Navigation across 4 billion miles** - hitting a 22-mile target after 13 years
5. **Scientific return** - revealing secrets of solar system formation

This project lets you **experience that discovery moment** firsthand - the thrill of finding something new in the vast darkness of space.

## 🔗 Learn More

- [NASA New Horizons Mission](https://www.nasa.gov/mission_pages/newhorizons/main/index.html)
- [Arrokoth Discovery Paper](https://science.nasa.gov/solar-system/kuiper-belt/arrokoth-2014-mu69/)
- [MAST Archive](https://archive.stsci.edu/) - Source of HST data
- [Hubble Space Telescope](https://hubblesite.org/)

---

_"The detailed images and other data that New Horizons could obtain from a KBO flyby will revolutionize our understanding of the Kuiper Belt and KBOs."_ - New Horizons Science Team

**Ready to make your own cosmic discovery? Fire up the scripts and find Arrokoth! 🌌**
