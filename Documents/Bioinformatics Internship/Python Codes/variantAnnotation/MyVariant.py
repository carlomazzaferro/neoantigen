
'''
Python Client for MyVariant.Info services
'''
from __future__ import print_function
import sys
import time
from itertools import islice
from collections import Iterable

import requests
try:
    from pandas import DataFrame
    from pandas.io.json import json_normalize
    df_avail = True
except:
    df_avail = False

try:
    import requests_cache
    caching_avail = True
except:
    caching_avail = False

__version__ = '0.2.0'

if sys.version_info[0] == 3:
    str_types = str
else:
    str_types = (str, unicode)



class MyVariantInfo:
    '''This is the client for MyVariant.info web services.
    Example:
        >>> mv = MyVariantInfo()
    '''
    def __init__(self, url='http://myvariant.info/v1'):
        self.url = url
        if self.url[-1] == '/':
            self.url = self.url[:-1]
        self.max_query = 1000
        # delay and step attributes are for batch queries.
        self.delay = 1
        self.step = 1000
        self.scroll_size = 1000
        # raise requests.exceptions.HTTPError for status_code > 400
        #   but not for 404 on getvariant
        #   set to False to surpress the exceptions.
        self.raise_for_status = True
        self._cached = False


    def _get(self, url, params={}, none_on_404=False):
        debug = params.pop('debug', False)
        return_raw = params.pop('return_raw', False)
        headers = {'user-agent': "Python-requests_myvariant.py/%s (gzip)" % requests.__version__}
        res = requests.get(url, params=params, headers=headers)
        if debug:
            return res
        if none_on_404 and res.status_code == 404:
            return None
        if self.raise_for_status:
            # raise requests.exceptions.HTTPError if not 200
            res.raise_for_status()
        if return_raw:
            return res.text
        ret = res.json()
        if 'from_cache' in vars(res):
            ret['_from_cache'] = vars(res).get('from_cache')
        return ret

    def _post(self, url, params):
        return_raw = params.pop('return_raw', False)
        headers = {'content-type': 'application/x-www-form-urlencoded',
                   'user-agent': "Python-requests_myvariant.py/%s (gzip)" % requests.__version__}
        res = requests.post(url, data=params, headers=headers)
        if self.raise_for_status:
            # raise requests.exceptions.HTTPError if not 200
            res.raise_for_status()
        if return_raw:
            return res
        ret = res.json()
        if 'from_cache' in vars(res):
            ret['_from_cache'] = vars(res).get('from_cache')
        return ret

    def _format_list(self, a_list, sep=','):
        if isinstance(a_list, (list, tuple)):
            _out = sep.join([safe_str(x) for x in a_list])
        else:
            _out = a_list     # a_list is already a comma separated string
        return _out



    @property
    def metadata(self):
        '''Return a dictionary of MyVariant.info metadata.
        Example:
        >>> metadata = mv.metadata
        '''
        _url = self.url+'/metadata'
        return self._get(_url)

    def set_caching(self, cache_db='myvariant_cache', **kwargs):
        ''' Installs a local cache for all requests.
            **cache_db** is the path to the local sqlite cache database.'''
        if caching_avail:
            requests_cache.install_cache(cache_name=cache_db, **kwargs)
            self._cached = True
        else:
            print("Error: The requests_cache python module is required to use request caching.")
            print("See - https://requests-cache.readthedocs.io/en/latest/user_guide.html#installation")
        return

    def stop_caching(self):
        ''' Stop caching.'''
        if self._cached:
            requests_cache.uninstall_cache()
            self._cached = False
        return

    def clear_cache(self):
        ''' Clear the globally installed cache. '''
        try:
            requests_cache.clear()
        except:
            pass


    def _getvariants_inner(self, vids, **kwargs):
        _kwargs = {'ids': self._format_list(vids)}
        _kwargs.update(kwargs)
        _url = self.url + '/variant/'
        return self._post(_url, _kwargs)

    def getvariants(self, vids, fields=None, **kwargs):
        '''Return the list of variant annotation objects for the given list of hgvs-base varaint ids.
        This is a wrapper for POST query of "/variant" service.
        :param ids: a list/tuple/iterable or a string of comma-separated HGVS ids.
                    `More about hgvs id <http://docs.myvariant.info/en/latest/doc/data.html#id-field>`_.
        :param fields: fields to return, a list or a comma-separated string.
                       If not provided or **fields="all"**, all available fields
                       are returned. See `here <http://docs.myvariant.info/en/latest/doc/data.html#available-fields>`_
                       for all available fields.
        :param as_dataframe: if True or 1 or 2, return object as DataFrame (requires Pandas).
                                  True or 1: using json_normalize
                                  2        : using DataFrame.from_dict
                                  otherwise: return original json
        :param df_index: if True (default), index returned DataFrame by 'query',
                         otherwise, index by number. Only applicable if as_dataframe=True.
        :return: a list of variant objects or a pandas DataFrame object (when **as_dataframe** is True)
        :ref: http://docs.myvariant.info/en/latest/doc/variant_annotation_service.html.

        .. Hint:: A large list of more than 1000 input ids will be sent to the backend
                  web service in batches (1000 at a time), and then the results will be
                  concatenated together. So, from the user-end, it's exactly the same as
                  passing a shorter list. You don't need to worry about saturating our
                  backend servers.
        .. Hint:: If you need to pass a very large list of input ids, you can pass a generator
                  instead of a full list, which is more memory efficient.
        '''
        if isinstance(vids, str_types):
            vids = vids.split(',') if vids else []
        if (not (isinstance(vids, (list, tuple, Iterable)))):
            raise ValueError('input "vids" must be a list, tuple or iterable.')
        if fields:
            kwargs['fields'] = self._format_list(fields)
        verbose = kwargs.pop('verbose', True)
        dataframe = kwargs.pop('as_dataframe', None)
        df_index = kwargs.pop('df_index', True)
        if dataframe in [True, 1]:
            dataframe = 1
        elif dataframe != 2:
            dataframe = None
        return_raw = kwargs.get('return_raw', False)
        if return_raw:
            dataframe = None

        query_fn = lambda vids: self._getvariants_inner(vids, **kwargs)
        out = []
        for hits in self._repeated_query(query_fn, vids, verbose=verbose):
            if return_raw:
                out.append(hits)   # hits is the raw response text
            else:
                out.extend(hits)
        if return_raw and len(out) == 1:
            out = out[0]
        if dataframe:
            out = self._dataframe(out, dataframe, df_index=df_index)
        return out


    def _fetch_all(self, **kwargs):
        ''' Function that returns a generator to results.  Assumes that 'q' is in kwargs.'''
        # get the total number of hits and start the scroll_id
        _url = self.url + '/query'
        res = self._get(_url, kwargs)
        try:
            scroll_id = res['_scroll_id']
            total_hits = int(res['total'])
        except KeyError:
            raise ScanError("Unable to open scroll.")
        if total_hits == 0:
            raise StopIteration
        kwargs.pop('q', None)
        kwargs.pop('fetch_all', None)
        print("Fetching {} variant(s)...".format(total_hits))
        while True:
            for hit in res['hits']:
                yield hit
            # get next scroll results
            kwargs.update({'scroll_id': scroll_id})
            res = self._get(_url, kwargs)
            if 'error' in res:
                break
            if '_warning' in res:
                print(res['_warning'])
            scroll_id = res.get('_scroll_id')

