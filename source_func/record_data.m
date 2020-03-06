function data = record_data(U, sensor_start, sensor_end)
x_start = sensor_start(1);
y_start = sensor_start(2);
x_end = sensor_end(1);
y_end = sensor_end(2);
data = [U(x_start,y_start:y_end-1)';U(x_start:x_end-1,y_end);...
    U(x_end,y_end:-1:y_start+1)';U(x_end:-1:x_start+1,y_start)];
end
