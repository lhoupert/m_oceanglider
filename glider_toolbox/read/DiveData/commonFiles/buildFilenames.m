% ****************************************************************************************************
%
% ****************************************************************************************************
function [filenames] = buildFilenames(folder, seagliderID, diveNumber)

    filenames.root      = sprintf('%3.3d%4.4d', seagliderID, diveNumber);
    filenames.baseFile  = sprintf('p%3.3d%4.4d', seagliderID, diveNumber);
    filenames.p_log     = fullfile(folder, sprintf('p%d%4.4d.log',    seagliderID, diveNumber));
    filenames.p_eng     = fullfile(folder, sprintf('p%d%4.4d.eng',    seagliderID, diveNumber));
    filenames.ppc_a_eng = fullfile(folder, sprintf('ppc%d%4.4da.eng', seagliderID, diveNumber));
    filenames.ppc_b_eng = fullfile(folder, sprintf('ppc%d%4.4db.eng', seagliderID, diveNumber));
    filenames.nc        = fullfile(folder, sprintf('p%d%4.4d.nc',     seagliderID, diveNumber));
    filenames.nc_gz     = fullfile(folder, sprintf('p%d%4.4d.nc.gz',  seagliderID, diveNumber));
    filenames.prs_a_eng = fullfile(folder, sprintf('prs%d%4.4da.eng', seagliderID, diveNumber));
    filenames.prs_b_eng = fullfile(folder, sprintf('prs%d%4.4db.eng', seagliderID, diveNumber));
    filenames.pcp_a_prf = fullfile(folder, sprintf('pcp%d%4.4da.prf', seagliderID, diveNumber));
    filenames.pcp_b_prf = fullfile(folder, sprintf('pcp%d%4.4db.prf', seagliderID, diveNumber));
    filenames.sound     = fullfile(folder, sprintf('p%d%4.4d_sound.csv', seagliderID, diveNumber));
end
