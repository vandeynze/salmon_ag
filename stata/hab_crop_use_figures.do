use "C:\Users\Lisa.Pfeiffer\Documents\GitHub\salmon_ag\output\land_use_habitat_use.dta", clear
drop if hab_use==.
reshape wide  count freq , i(name land_use Crop) j(hab_use)
order name land_use Crop cropgroup count* freq* 
recode count* freq* (.=0)
collapse (sum) count* freq*, by(name cropgroup)
egen total_freq=rsum(freq1 freq2 freq3)
drop if total_freq<0.01
drop if cropgroup=="other"|cropgroup=="Water"
gsort name -freq1
format freq* total_freq %9.3f
drop count*


gen species="Chinook" if strpos(name, "Chinook")
replace species="Coho" if strpos(name, "Coho")
replace species="Sockeye" if strpos(name, "Sockeye")
replace species="Pink" if strpos(name, "Pink")
replace species="Chum" if strpos(name, "Chum")
gen state="CA"  if strpos(name, "CA")|strpos(name, "California")|strpos(name, "Central Valley")|strpos(name, "Sacramento")
replace state="OR"  if strpos(name, "Columbia")|strpos(name, "Willamette")
graph bar (asis) freq1 freq2 freq3 if species=="Chinook" & state=="CA", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Chinook" & state=="OR", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Chinook" & state=="", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Coho", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Sockeye", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Pink", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
graph bar (asis) freq1 freq2 freq3 if species=="Chum", over(cropgroup) stack by(name) horizontal legend(order(1 "Spawning" 2 "Rearing" 3 "Migration") rows(1))
