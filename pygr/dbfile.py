
import shelve, anydbm, sys, UserDict
import logger

class WrongFormatError(IOError):
    'attempted to open db with the wrong format e.g. btree vs. hash'
    pass
class NoSuchFileError(IOError):
    'file does not exist!'
    pass
class PermissionsError(IOError):
    'inadequate permissions for requested file'
    pass
class ReadOnlyError(PermissionsError):
    'attempted to open a file for writing, but no write permission'
    pass

def open_anydbm(*args, **kwargs):
    'trap anydbm.error message and transform to our consistent exception types'
    try:
        return anydbm.open(*args, **kwargs)
    except anydbm.error, e:
        msg = str(e)
        if msg.endswith('new db'):
            raise NoSuchFileError(msg)
        elif msg.startswith('db type'):
            raise WrongFormatError(msg)
        raise

try: # detect whether bsddb module available and working...
    import bsddb
    try:
        bsddb.db
    except AttributeError:
        raise ImportError
except ImportError:
    bsddb = None

def open_bsddb(filename, flag='r', useHash=False, mode=0666):
    """open bsddb index instead of hash by default.
    useHash=True forces it to use anydbm default (i.e. hash) instead.
    Also gives more meaningful error messages."""
    try: # 1ST OPEN AS BTREE
        if useHash: # FORCE IT TO USE HASH INSTEAD OF BTREE
            return open_anydbm(filename, flag)
        else:
            return bsddb.btopen(filename, flag, mode)
    except bsddb.db.DBAccessError: # HMM, BLOCKED BY PERMISSIONS
        if flag=='c' or flag=='w': # TRY OPENING READ-ONLY
            try:
                ifile = file(filename)
            except IOError: # HMM, NOT EVEN READABLE. RAISE GENERIC PERMISSIONS ERROR
                raise PermissionsError('insufficient permissions to open file: '
                                       +filename)
            ifile.close() # OK, WE CAN READ THE FILE, SO RAISE EXCEPTION WITH
            raise ReadOnlyError('file is read-only: '+filename) # VERY SPECIFIC MEANING!
        else: # r OR n FLAG: JUST RAISE EXCEPTION
            raise PermissionsError('insufficient permissions to open file: '
                                   +filename)
    except bsddb.db.DBNoSuchFileError:
        raise NoSuchFileError('no file named: '+filename)
    except bsddb.db.DBInvalidArgError: # NOT A BTREE FILE...
        try:
            if useHash: # NO POINT IN TRYING HASH YET AGAIN...
                raise bsddb.db.DBInvalidArgError
            # fallback to using default: hash file
            return open_anydbm(filename, flag)
        except bsddb.db.DBInvalidArgError:
            raise WrongFormatError('file does not match expected shelve format: '+filename)

def open_index(filename, flag='r', useHash=False, mode=0666):
    if bsddb is None:
        d = open_anydbm(filename, flag)
        if not useHash:
            logger.warn('Falling back to hash index: unable to import bsddb')
        return d
    return open_bsddb(filename, flag, useHash, mode)

def iter_gdbm(db):
    'iterator for gdbm objects'
    k = db.firstkey()
    while k is not None:
        yield k
        k = db.nextkey(k)

class _ClosedDict(UserDict.DictMixin):
    """This dummy class exists solely to raise a clear error msg if accessed.
    Copied from the Python 2.6 shelve.py """
    def closed(self, *args):
        raise ValueError('invalid operation on closed shelf')
    __getitem__ = __setitem__ = __delitem__ = keys = closed

    def __repr__(self):
        return '<Closed Dictionary>'



class BetterShelf(shelve.Shelf):
    """Shelf subclass that fixes its horrible iter implementation.
    """
    def __iter__(self):
        'avoid using iter provided by shelve/DictMixin, which loads all keys!'
        try:
            return iter(self.dict)
        except TypeError: # gdbm has wierd iterator behavior...
            return iter_gdbm(self.dict)

    if sys.version_info < (2, 6): # Python finally added good err msg in 2.6
        def close(self):
            if isinstance(self.dict, _ClosedDict):
                return # if already closed, nothing further to do...
            shelve.Shelf.close(self) # close Shelf as usual
            self.dict = _ClosedDict() # raises sensible error msg if accessed

def shelve_open(filename, flag='c', protocol=None, writeback=False,
                useHash=True, mode=0666, *args, **kwargs):
    """improved implementation of shelve.open() that won't generate
bogus __del__ warning messages like Python's version does."""
    d = open_index(filename, flag, useHash, mode) # construct Shelf only if OK
    return BetterShelf(d, protocol, writeback, *args, **kwargs)
