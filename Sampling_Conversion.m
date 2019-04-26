abs(mean(AOUT_Orig_Data(134318:134318+5e3) - AOUT_Trojan_Data(135862:135862+5e3)))

std(diff(AOUT_Orig_Data(134318:134318+5e3) - AOUT_Trojan_Data(135862:135862+5e3)))


abs(mean(AOUT_Orig_Data(134318:134318+5e3) - AOUT_Defense_Data(133182:133182+5e3)))

std(diff(AOUT_Orig_Data(134318:134318+5e3) - AOUT_Defense_Data(133182:133182+5e3)))

% Creating Constant Time Steps

AOUT_Orig_Time = linspace(min(AOUT_Orig_Time),max(AOUT_Orig_Time),length(AOUT_Orig_Time));

AOUT_Orig_Time = AOUT_Orig_Time';

AOUT_Trojan_Time = linspace(min(AOUT_Trojan_Time),max(AOUT_Trojan_Time),length(AOUT_Trojan_Time));

AOUT_Trojan_Time = AOUT_Trojan_Time';

AOUT_Defense_Time = linspace(min(AOUT_Defense_Time),max(AOUT_Defense_Time),length(AOUT_Defense_Time));

AOUT_Defense_Time = AOUT_Defense_Time';

% Performing Conversion

[P1,Q1] = rat((1/(AOUT_Orig_Time(2)-AOUT_Orig_Time(1)))/(1/(AOUT_Defense_Time(2)-AOUT_Defense_Time(1))));

AOUT_Defense_Data = resample(AOUT_Defense_Data,P1,Q1);

[P2,Q2] = rat((1/(AOUT_Orig_Time(2)-AOUT_Orig_Time(1)))/(1/(AOUT_Trojan_Time(2)-AOUT_Trojan_Time(1))));

AOUT_Trojan_Data = resample(AOUT_Trojan_Data,P2,Q2);
