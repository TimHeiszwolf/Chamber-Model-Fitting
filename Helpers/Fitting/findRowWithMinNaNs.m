function rowIndex = findRowWithMinNaNs(inputMatrix)
    % findRowWithMinNaNs Finds the index of the row with the least number of NaNs.
    % This function is pretty useless/uncececarry. It is only used for when importing Thomson scattering (which is not really used). Kept it here for if ever needed but can safely disregards.    
    nanCountsPerRow = sum(isnan(inputMatrix), 2);

    [~, rowIndex] = min(nanCountsPerRow);
end