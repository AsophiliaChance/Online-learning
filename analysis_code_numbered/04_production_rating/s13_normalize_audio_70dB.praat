form Batch normalize wav files to 70 dB
    sentence input_folder   D:\online_learning\production\syllable\norm70db\test\
    sentence output_folder  D:\online_learning\production\syllable\norm70db\test2\
    real target_dB 70
endform

# ---------- normalize folder path ----------
if right$(input_folder$, 1) <> "\" and right$(input_folder$, 1) <> "/"
    input_folder$ = input_folder$ + "\"
endif
if right$(output_folder$, 1) <> "\" and right$(output_folder$, 1) <> "/"
    output_folder$ = output_folder$ + "\"
endif

# ---------- list wav files ----------
Create Strings as file list: "wavList", input_folder$ + "*.wav"
num_files = Get number of strings

if num_files = 0
    printline "No wav files found in: " + input_folder$
    exit
endif

# ---------- main loop ----------
for i from 1 to num_files
    select Strings wavList
    wav_name$ = Get string: i
    full_path$ = input_folder$ + wav_name$

    # Read sound
    Read from file: full_path$
    sound_id = selected("Sound")

    # Normalize and save
    Scale intensity: target_dB
    Save as WAV file: output_folder$ + wav_name$

    # Remove sound safely (by id)
    selectObject: sound_id
    Remove
endfor

printline ""
printline "Done."
printline "Processed files: " + string$(num_files)
printline "Target intensity: " + string$(target_dB) + " dB"
