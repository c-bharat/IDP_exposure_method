
/* 

CODE TO FACILITATE IMPLEMENTATION OF THE INDIVIDUALISED DISPENSING PATTERNS (IDP) METHOD FOR DEFINING MEDICINE EXPOSURE

Drafted by: Chrianna Bharat with input from Luke Buizen
Ref citation:  Bharat, C, Degenhardt, L, Pearson, S-A, et al. A data-informed approach using individualised dispensing patterns to estimate medicine exposure periods and dose from pharmaceutical claims data. Pharmacoepidemiol Drug Saf. 2023; 32( 3): 352- 365. doi:10.1002/pds.5567

INPUT DATA: 
individual-level daily dispensing data, where dispensings on the same day for the same medicine(s) are combined/summed together (daily total dispensed)
(called macro_d in the code below)

NECESSARY VARIABLES in input data: 
- PPN: person-level ID
- group: a variable identifying all dispensings for a given medicine/medicine class 
- date_of_supply: date of dispensing d_n
- q_D: the quantity dispensed at d_n (this may be recorded differently in different datasets e.g., number of pills v number of box/s - be sure to consider the appropriateness and whether any modifications are needed prior to running the below)
- Date_of_Supply_index: date of first d_n
- DateDeath: date of death

Important: in the following code the study end of follow-up/right censoring date was '31DEC2018'D - this should be updated in each project 
	
OUTPUT: 
change history data where each individual's follow-up time is divided into current, recent and formally exposed intervals
(this dataset is named based on the 'output' macro variable in the exposure_by_drug MACRO function)

Key variables include...
- start_date: 1-date of the start of the interval 
- end_date: date of the end of the interval 
- cens_date: censoring date (e.g., date of death, end of follow-up) 
- es: 1 = currently exposed, 2 = recently exposed, 3 = formerly exposed (to the given medicine/medicine class)

*/

* STEP 1: Define population-based exposure estimate; 
* Needs to be run for each medicine/medicine class separately - called by &item_code; 
%MACRO e_pop_estimate(item_code);

	PROC SORT DATA = macro_d; 
		BY PPN Date_of_Supply; 
	RUN; 

	DATA e_pope new_episode;
		SET macro_d(WHERE = (group=&item_code.));
		RETAIN last_time last_q t_nm1 q_nm1 new_episode;
		BY PPN Date_of_Supply;

		IF first.PPN THEN DO;
			t_nm1=.;
			q_nm1=.;
			last_time=.;
			last_q=.;
			new_episode=0;
		END;

		t_nm1=last_time;
		q_nm1=last_q;
		last_q=q_D;
		last_time=Date_of_Supply;
		IF t_nm1=. | Date_of_Supply-t_nm1>=365 THEN new_episode=1;
											   ELSE OUTPUT e_pope;
		IF new_episode=1 THEN DO;
			OUTPUT new_episode;
			new_episode=0;
		END;
		FORMAT last_time date9. t_nm1 date9.;
	RUN;
	PROC PRINT DATA = macro_d (OBS=50); WHERE group=&item_code.; RUN; 
	PROC PRINT DATA = e_pope (OBS=50); RUN; 
	PROC PRINT DATA = new_episode (OBS=50); RUN; 

	*Calculate the number of days per unit following a dispensing day;
	DATA e_pop_est;
		SET e_pope;
		dif = Date_of_Supply-t_nm1; 
		e_pop_est = dif/q_nm1;

		BY ppn; 
		IF first.ppn then first_interval = 1; 
	RUN;
	PROC PRINT DATA = e_pop_est (OBS=50); RUN; 

	* Calculate the 80th percentile estimate of the distribution for days per unit;
	PROC UNIVARIATE DATA = e_pop_est ; 
		VAR e_pop_est; 
		OUTPUT 	OUT = E_est_item_code
				N =N
				PCTLPRE = P_ 
				PCTLPTS=20,50,80,90; 
	RUN;
	DATA E_est_item_code; 
		SET E_est_item_code; 
		item_code="&item_code."; 
	RUN;

	DATA combined_item_code;
		SET combined_item_code
			E_est_item_code;
	RUN;
	PROC PRINT DATA = E_est_item_code (OBS=50); RUN; 
	PROC PRINT DATA = combined_item_code (OBS=50); RUN; 
	
%MEND e_pop_estimate;

* Empty tables for population-level days per quantity per day; 
DATA combined_item_code;
	LENGTH item_code $6. N 8. P_20 8. P_50 8. P_80 8. P_90 8.;
RUN;

/* Run e_pop_estimate macro */ 

* Days per quantity estimates; 
DATA Combined_item_code3;
	SET Combined_item_code;
	LABEL p_20 = "20th percentile";
	LABEL p_50 = "50th percentile";
	LABEL p_80 = "80th percentile";
	LABEL p_90 = "90th percentile";
	LABEL N = "N";

	group  = input(item_code,8.);
RUN;

/* Create change history file using the IDP method */ 
* Note that the macro variable drug_number is a drug group reference number (project_specific and potentially useful when exposure to multiple different medicines are being considered); 
* Needs to be run for each medicine/medicine class separately - called by &drug_number; 
%macro exposure_by_drug(drug_number, output);

	* Subset dispensing-day data to medicines/dispensings of interest;
	DATA macro_d_temp; 
		SET macro_d; 
		IF group=&drug_number. AND date_of_supply>=Date_of_Supply_index;  
		all=1; 
	RUN;
	PROC PRINT DATA = macro_d (OBS=5); RUN;
	PROC PRINT DATA = macro_d_temp (OBS=5); RUN; 
	
	* Merge in population estimate; 
	DATA E_est_macro; 
		SET Combined_item_code3(WHERE=(group=&drug_number.) 
									KEEP=group P_80); 
		all=1; 
	RUN;
	PROC PRINT DATA = E_est_macro (OBS=5); RUN;
	DATA macro_d_temp; 
		MERGE macro_d_temp 
			  E_est_macro; 
		BY all; 
	RUN;

	* Implement IDP evaluation of duration of exposure; 
	DATA macro_d1;
		SET macro_d_temp ;
		RETAIN 	ep_num episode_dispensing recent_exp 
				t_nm1 t_nm2 t_nm3 q_nm1 q_nm2 q_nm3 
				e_n 
				last_time last_q 
				first_&drug_number._date
				unique_ep_id;
		FORMAT last_time DATE9.;

		BY PPN Date_of_Supply;

		recent_exp=7; * Specify the max length of recent exposure periods here; 

		* For the first obs, set all parameters to zero/missing; 
		IF first.PPN THEN DO;
			ep_num=0;
			episode_dispensing=0;
			t_nm1=.;
			t_nm2=.;
			t_nm3=.;
			q_nm1=.;
			q_nm2=.;
			q_nm3=.;
			e_n = 99999999;
			first_&drug_number._date=date_of_supply; 
		END;
		no_formerly_exposed=0;
		
		* Set paramaters to contingent on it being >= nth dispensing for that episode; 
		IF episode_dispensing>=3 THEN t_nm3=t_nm2;
		IF episode_dispensing>=2 THEN t_nm2=t_nm1;
		IF episode_dispensing>=1 THEN t_nm1=last_time;

		IF episode_dispensing>=3 THEN q_nm3=q_nm2;
		IF episode_dispensing>=2 THEN q_nm2=q_nm1;
		IF episode_dispensing>=1 THEN q_nm1=last_q;

		* This looks to see if the gap between dispensings is short enought to be from the same episode;
		IF ep_num = 0 OR (Date_of_Supply > (t_nm1 + e_n + recent_exp)) THEN DO;
			ep_num+1;
			unique_ep_id+1;
			episode_dispensing=1;
			t_nm1=.;
			t_nm2=.;
			t_nm3=.;
			q_nm1=.;
			q_nm2=.;
			q_nm3=.;
			e_n = q_D*(P_80);
		END;
		ELSE DO; 
			* If the interval contains at least one day while formerly exposed; 
			IF Date_of_Supply > (t_nm1 + e_n + recent_exp) THEN DO;
				ep_num+1; 
				no_formerly_exposed=1; 
			END;
			episode_dispensing+1;

			* STEP 2 / 3: Calculate the estimated number of days of exposure; 
				 IF t_nm1=. THEN e_n = q_D*(P_80);
			ELSE IF t_nm2=. THEN e_n = q_D*( (3/6)*(Date_of_Supply-t_nm1)/q_nm1 + (2/6)*(P_80) 				+ (1/6)*(P_80) 			);
			ELSE IF t_nm3=. THEN e_n = q_D*( (3/6)*(Date_of_Supply-t_nm1)/q_nm1 + (2/6)*(t_nm1-t_nm2)/q_nm2 + (1/6)*(P_80) 			);
			ELSE 				 e_n = q_D*( (3/6)*(Date_of_Supply-t_nm1)/q_nm1 + (2/6)*(t_nm1-t_nm2)/q_nm2 + (1/6)*(t_nm2-t_nm3)/q_nm3	);

		END;

		last_q=q_D;
		last_time=Date_of_Supply;

	RUN;
	PROC PRINT DATA = macro_d1 (OBS=50); RUN; 

	* Merge death dates with PBS data; 
	DATA macro_d1;
		FORMAT cens_date DATE9.; 
		MERGE macro_d1(IN=keep) 
			  clean.ndi(KEEP=PPN DeathDate 
						IN=inndi);
		BY PPN;
		IF inndi=1 THEN death=1; 
				   ELSE death=0;
		IF keep=1;
		cens_date = MIN(DeathDate, '31DEC2018'D);
	RUN;
	
	* Look ahead to get date of next dispensing; 
	PROC SORT DATA=macro_d1; 
		BY PPN DESCENDING Date_of_Supply; 
	RUN;
	DATA macro_d2;
		SET macro_d1;

		PPN1=lag(PPN);
		SEE1=lag(Date_of_Supply);
		EP1=lag(ep_num);
		FORMAT SEE1 date9.;

		BY ppn; 
		IF first.ppn THEN DO; 
			SEE1 = '31DEC9999'D; * Used as an extreme end point - this creates a row from the censoring date to the extreme data, and should be dropped before analysis; 
			last = 1; 
			EP1 = .; 
		END; 
	RUN;

	* Define intervals of current, recent and former exposure - includes Step 4; 
	PROC SORT DATA=macro_d2; 
		BY PPN Date_of_Supply; 
	RUN; 
	DATA macro_episodes
		 post_death; 
		RETAIN start_date end_date; 
		FORMAT start_date end_date DATE9.;
		SET macro_d2; 

		* Currently exposed; 
		es=1; 	
		start_date = date_of_supply-1; 
		end_date = MIN(INTNX('DAY',date_of_supply,e_n,'END'), SEE1-1, DeathDate, '31DEC2018'D);

		* If post_death then delete ; 
		IF date_of_supply > MIN(DeathDate, '31DEC2018'D) THEN DO; 
			OUTPUT post_death; 
		END; 
		ELSE DO; 
			OUTPUT macro_episodes; 

			*Recently exposed; 
			IF MIN(SEE1-1,DeathDate,'31DEC2018'D) > INTNX('DAY',date_of_supply,e_n,'END') THEN DO; 
				es=2;
				start_date = INTNX('DAY',date_of_supply,e_n,'END');  
				end_date = MIN(INTNX('DAY',date_of_supply,e_n+recent_exp,'END'), SEE1-1, DeathDate, '31DEC2018'D);
				OUTPUT macro_episodes; 

				* Formerly exposed; 
				IF MIN(SEE1-1,DeathDate,'31DEC2018'D) > INTNX('DAY',date_of_supply,e_n+recent_exp,'END') THEN DO; 
					es=3;
					start_date = INTNX('DAY',date_of_supply,e_n+recent_exp,'END'); 
					end_date = MIN(SEE1-1, DeathDate, '31DEC2018'D);
					OUTPUT macro_episodes; 
				END; 
			END; 
		END; 
	RUN; 

	* Clean up variables and specify those wanted for final output; 
	PROC SORT DATA = macro_episodes; BY PPN start_date; RUN; 

	* Create intervals based on days since start of follow-up and carry episode start date for every subsequent row until start of next episode; 
	DATA macro_episodes2; 
		RETAIN ep_st_date ; 
		FORMAT ep_st_date DATE9.;
		SET macro_episodes;  

		IF es = 2 or es = 3 THEN DO; 
			date_of_supply = ''D; 
			n_D=.;
			q_D=.; 
			e_n=.;
			episode_dispensing = .; 
		END; 
		
		start_t = start_date - (Date_of_Supply_index365-1); 
		end_t = end_date - (Date_of_Supply_index365-1); 
		
		BY ppn ep_num;
		IF first.ep_num THEN ep_st_date = date_of_supply; 
	RUN; 					

	* Create episode number for each persons exposure interval; 
	DATA &output.; 
		RETAIN ep_num&drug_number. first_&drug_number._date; 
		FORMAT first_&drug_number._date DATE9.;
		SET macro_episodes2 (DROP=ep_num);   

		pdays = end_date - start_date;
		lag_es = LAG(es);
		
		BY ppn ;
		IF first.ppn THEN DO; 
			ep_num&drug_number. = 1; 
			rec_num=1; 
			first_&drug_number._date = Date_of_Supply; 
			lag_es=.; 
		END; 
		ELSE IF lag_es NE es THEN rec_num+1; 

		* New episode defined following at least one day of 'formerly exposed'; 
		IF es=1 AND lag_es=3 THEN ep_num&drug_number.+1; 
	RUN; 		

%MEND exposure_by_drug;

* Every dataset is different - be sure to conduct logic checks on output dataset prior to using in any analyses; 
