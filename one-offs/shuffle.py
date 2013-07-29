def shuffle_leptons(event):
	# make an OF pair if possible
	if abs(event['pdg1']) == abs(event['pdg2']):
		if abs(event['pdg1']) != abs(event['pdg3']):
			event['pdg2'], event['pdg3'] = event['pdg3'], event['pdg2']
			event['pt2'], event['pt3'] = event['pt3'], event['pt2']
			event['eta2'], event['eta3'] = event['eta3'], event['eta2']
			event['phi2'], event['phi3'] = event['phi3'], event['phi2']

		elif abs(event['pdg2']) != abs(event['pdg3']):
			event['pdg1'], event['pdg3'] = event['pdg3'], event['pdg1']
			event['pt1'], event['pt3'] = event['pt3'], event['pt1']
			event['eta1'], event['eta3'] = event['eta3'], event['eta1']
			event['phi1'], event['phi3'] = event['phi3'], event['phi1']

	return event