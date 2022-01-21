
primaryFromISETBioGaborCal = screenCalObjFromISETBio.get('P_device') \ ...
    (ISETBioGaborCal-screenCalObjFromISETBio.get('P_ambient'));
settingsFromISETBioGaborCal = PrimaryToSettings(screenCalObjFromISETBio,primaryFromISETBioGaborCal);
if (max(abs(standardSettingsGaborCal(:)-settingsFromISETBioGaborCal(:))./standardSettingsGaborCal(:)) > 1e-6)
    error('Cannot get home again in settings land');
end