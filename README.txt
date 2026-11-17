[This README.txt template was created by Bob Montgomerie (Queen’s University) and was  updated on 08 Oct 2024; please replace or delete instructions in square brackets]


GENERAL INFORMATION

1. Paper Citation
Do sperm whalee responses during interspecies interactions reflect potential threats? (submitted)

2. Brief abstract

Predation pressure is a key driver of sociality, as aggregating confers both passive and active defense against predators. For sperm whales (Physeter macrocephalus), a deep-diving and nomadic species, cooperative protection of vulnerable young against predators may be a primary benefit of maintaining social bonds. Here, we determine whether sperm whale collective behaviours during interspecific interactions correspond with threat levels and social context. We examined aerial footage of sperm whale groups collected off the Galapagos Islands in 2023 before, during and after interactions with different odontocete species, including Tursiops spp., Orcinus orca, Delphinus spp., and Globicephala macrorhynchus (n=12 interactions, 11 days). We analyzed sperm whales’ spatial configuration and cluster size across contexts and assessed whether the relative position of calves shifts in response to varying risk levels. We report the first case of sperm whales collectively chasing away killer whales, providing empirical support for the hypothesis that variation in sperm whale group size across ocean basins is driven by predation pressure. Additionally, we show that sperm whales bunched more tightly and with more variable headings during aversive interactions. Although calves were more likely to be found towards the centre-front of clusters, their relative positions did not shift during interactions. By analyzing the collective response of sperm whales to varying levels of threat, we provide evidence for the adaptive nature of the species' anti-predator behaviour and highlight the role that cooperative defense plays in sperm whale society. 

3. Originators

Ana Eguiguren*, Luke Rendell, Elizabeth Zwamborn, Isabel Fernandez-McAuley, Hal Whitehead

4. Contact information

*anaeguibur@gmail.com

5. Date of data collection

February - May, 2023

6. Geographic location(s) of data collection

Galapagos Marine Reserve and surrounding waters

7. Information about funding sources that supported the collection and curation of the
data

Data was collected with support from Canada's National Science and Engineeing Research Council (NSERC) and
the Rufford Foundation. 

ACCESS INFORMATION

1. Licenses/restrictions placed on the data

CC-BY-NC-SA 4.0 Attribution-NonCommercial-ShareAlike


DATA FILES AND VARIABLES


1. labled_instances.csv
Each row contains the positions of an individual whale within a frame.

 • date: YYYY-MM-DD
 • group_ID: sperm whale group identity assign following Whitehead & Whitehead (1997)
 • frame_local_time: frame timestamp in local time (UTC-6) YYYY-MM-DDD HH:MM:SS.ms
 • interaction: interaction number (identifies point in time when other species was withing <30 m from a sperm whale)
 • event: event ID (identified as a set of consecutive recordings where other species were observed interacting with sperm whales)
 • video: full path to video recording with original file name
 • video_file: file name of video recording
 • frame_idx: frame number within recording
 • instance_id: identifies individual whales labelled within frames
 • frame_id: uniquely identifies frames accross recordings
 • snapshot path: location of stored snapshot
 • x_snout: horizontal position of sperm whale's snout
 • y_snout: vertical position of sperm whale's snout
 • x_dorsal: horizontal position of sperm whale's dorsal fin
 • y_dorsal: vertical position of sperm whale's dorsal fin
 • x_fluke: horizontal position of sperm whale's fluke notch
 • y_fluke: vertical position of sperm whale's fluke notch
 • age_class: coarse age classification of individuals based on 'good measurements'
 • Intensity: aversity level (Neutral vs. Aversive)
 • species: species of other interacting individuals (KW = killer whales; Gm = short-finned pilot whales; Tt = Tursiops truncatus/Delphinus spp.)







2. summary_5min_behavioural_state.csv

[brief description of the file followed by a bulleted listing of each variable/column with a detailed description]

   • variable: description
   • variable: description
   • variable: description
 

CODE SCRIPTS AND WORKFLOW

[list each code script file with a brief description; if there is more than one code script, 
provide information about the best workflow for running those files]

1. filename

[brief description]

2. filename

[brief description]


SOFTWARE VERSIONS

[list all software used with version numbers of the software and any loaded packages]

software name and version

loaded packages:

 • package: version
 • package: version
 • package: version
 • etc

software name and version

loaded packages:

 • package: version
 • package: version
 • package: version
 • etc

REFERENCES

[list all references mentioned in this README file]

