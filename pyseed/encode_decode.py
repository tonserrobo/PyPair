
class Encode_decode():
    '''Encoder/Decoder class for 2-bit encoding of chars'''
    def encode(self, index):
        '''encode the data set using 2bit encoding'''
        encoded =[]
        dictionary = {'$': '', ',':'', 'A': '00', 'C': '01', 'G': '10', 'T': '11'}

        if isinstance(index, list):
            # encode list of strings
            for i in range(len(index)):
                substring = index[i]
                transTable = substring.maketrans(dictionary)
                txt = substring.translate(transTable)
                encoded.append(txt)
        else:
            # encode single string
            transTable = index.maketrans(dictionary)
            txt = index.translate(transTable)
            encoded = txt
        return encoded

    def decode(self, index):
        '''Decode a string or list of strings from 2-bit encoding to their respective char'''
        decoded = []
        decoded_substring =[]
        dictionary = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
        
        if isinstance(index, list):
            # decode list of strings
            for i in range(len(index)):
                count = 0
                while count < len(index):
                    i = count
                    j = i + 1
                    substring = index[i] + index[j]
                    decoded_substring.append(dictionary.get(str(substring), str(substring)))
                    count +=2
                decoded_string = ''.join(decoded_substring)
                decoded.append(decoded_string)
        else:
            # decode single string
            count = 0
            while count < len(index):
                i = count
                j = i + 1
                substring = index[i] + index[j]
                decoded_substring.append(dictionary.get(str(substring), str(substring)))
                count +=2
            decoded_string = ''.join(decoded_substring)
            decoded = decoded_string
        return decoded