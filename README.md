# BahaviourAnalysis
Code to analise data from a 2AFC behaviour task

Each row corresponds to a separate trial.

Animal: Number assigned to each rat.
Session: Number of the session the trials corresponds to.
Block: Block index referenced to the current session.
Trial_start: Trial start time. (dd-mmm-yy hh-mm-ss)
Trial_end: Trial end time. (dd-mmm-yy hh-mm-ss)
Trial_duration Trial duration. (Seconds)
LED_frequency: LED frequency used as a stimulus in the current trial.

Success: 0: aborted trial. 1: correct trial -1: incorrect trial

timed_ITI: Time elapsed between the end of the previous trial and the start of the current trial. (Seconds)

LED_cue: Time the animal is forced to experience the LED stimulus. 1 second.

response_poke: 2: Left poke 3: Right poke NaN: no response

_______

Psychometric curve is calculated by fitting a sigmoid curve to the averaged responses per conditions. To calculate the standard deviation a bootstrap method was applied by resampling with replacment (50 iterations used). 

_______

Code for Pearson correlation of behavioural responses with a different data modality (in this case was fMRI) is avaialable.

