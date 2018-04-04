#!/usr/bin/env python

from matplotlib import pyplot as plt

radar_nis = []
lidar_nis = []
radar_threshold = 7.815;
lidar_threshold = 5.991;

with open('nis_data.csv', 'r') as fd:
	for line in fd:
		sensor_type, nis = [float(x) for x in line.split(',')]
		if int(sensor_type):
			radar_nis.append(nis)
		else:
			lidar_nis.append(nis)

ax = plt.subplot(211)
ax.axhline(radar_threshold, color="red");
plt.plot(radar_nis, "-")
plt.title("RADAR NIS")
plt.ylabel("NIS")
plt.xlabel("time")
plt.legend(["95% threshold", "RADAR NIS"])
plt.ylim((0, 15))
plt.xlim((0, len(radar_nis)))
ax = plt.subplot(212)
ax.axhline(lidar_threshold, color="red");
plt.plot(lidar_nis, "-")
plt.title("LIDAR NIS")
plt.ylabel("NIS")
plt.xlabel("time")
plt.legend(["95% threshold", "RADAR NIS"])
plt.ylim((0, 15))
plt.xlim((0, len(radar_nis)))
plt.subplots_adjust(hspace=0.59)
plt.savefig('NIS_plot.png')
