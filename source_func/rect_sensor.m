function sensor = rect_sensor(rect_start, rect_end, gridx, gridy)
sensor = [gridx(rect_start(1),rect_start(2):1:rect_end(2)-1);gridy(rect_start(1),rect_start(2):1:rect_end(2)-1)];
sensor = [sensor [gridx(rect_start(1):1:rect_end(1)-1,rect_end(2)) gridy(rect_start(1):1:rect_end(1)-1,rect_end(2))]'];
sensor = [sensor [gridx(rect_end(1),rect_end(2):-1:rect_start(1)+1);gridy(rect_end(1),rect_end(2):-1:rect_start(1)+1)]];
sensor = [sensor [gridx(rect_end(1):-1:rect_start(1)+1,rect_start(2)) gridy(rect_end(1):-1:rect_start(1)+1,rect_start(2))]'];
end