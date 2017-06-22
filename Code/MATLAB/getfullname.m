function fullName = getfullname(namePre, i, namePost, fileType)

if isnan(i)
    fullName = [namePre, namePost, '.', fileType];
else
    fullName = [namePre, num2str(i), namePost, '.', fileType];
end

end