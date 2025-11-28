num_files = 5; % Number of files to load
sonar_data = cell(num_files, 1); % Preallocate cell array to hold data

for i = 1:num_files
    filename = sprintf("Sonar_test_data_lyddodtrum/Sonar_test_data/Sonar%d_5m.dat", i);
    sonar_data{i} = load(filename); % Load each .dat file into the cell array
end
sonar0_5 = load("Sonar_test_data_lyddodtrum\Sonar_test_data\Sonar0_5m.dat");
plot(sonar0_5);
