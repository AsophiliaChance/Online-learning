form Batch cut by paired boundaries (no filters, no decimals in names)
    text input_dir  D:\online_learning\production\task_audio\
    text output_dir D:\online_learning\production\syllable\
    sentence tier_name  John
    real pad_s  0.00
endform

Create Strings as file list: "fileList", input_dir$ + "*.wav"
select Strings fileList
n = Get number of strings

for i from 1 to n
    select Strings fileList
    wav$ = Get string: i

    len = length (wav$)
    if len < 5
        goto nextFile
    endif
    base$ = left$ (wav$, len - 4)

    wavPath$ = input_dir$ + wav$
    tgPath$  = input_dir$ + base$ + ".TextGrid"

    if fileReadable (tgPath$) = 0
        goto nextFile
    endif

    Read from file: wavPath$
    Rename: "curSound"
    Read from file: tgPath$
    Rename: "curTG"

    select TextGrid curTG
    nt = Get number of tiers
    tier_index = 0
    for t from 1 to nt
        name$ = Get tier name: t
        if name$ = tier_name$
            tier_index = t
        endif
    endfor

    if tier_index = 0
        goto cleanupFile
    endif

    select TextGrid curTG
    nint = Get number of intervals: tier_index

    syl = 0
    nmax = nint - 1
    k = 1

    while k <= nmax
        select TextGrid curTG
        t1 = Get end time of interval: tier_index, k
        t2 = Get end time of interval: tier_index, k + 1

        if t2 > t1
            t1p = t1 - pad_s
            if t1p < 0
                t1p = 0
            endif
            t2p = t2 + pad_s

            select Sound curSound
            durAll = Get total duration
            if t2p > durAll
                t2p = durAll
            endif

            if t2p > t1p
                select Sound curSound
                Extract part: t1p, t2p, "rectangular", 1.0, "yes"
                Rename: "curPart"

                syl = syl + 1
                if syl <= 99
                    if syl < 10
                        num$ = "0" + string$(syl)
                    else
                        num$ = string$(syl)
                    endif

                    outPath$ = output_dir$ + base$ + "_" + num$ + ".wav"
                    select Sound curPart
                    Save as WAV file: outPath$
                endif

                select Sound curPart
                Remove
            endif
        endif

        k = k + 2
    endwhile

    printline "Saved: " + base$ + " = " + string$(syl)

label cleanupFile
    select TextGrid curTG
    Remove
    select Sound curSound
    Remove

label nextFile
endfor

select Strings fileList
Remove
printline "Done."
