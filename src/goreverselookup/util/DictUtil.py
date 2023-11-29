# Dictionary utility class

class DictUtil:
	@classmethod
	def reverse_dict(cls, input:dict):
		"""
		Reverses the keys and values in the input dictionary.
		"""
		res = {}
		for key,value in input.items():
			res[value] = key
		return res