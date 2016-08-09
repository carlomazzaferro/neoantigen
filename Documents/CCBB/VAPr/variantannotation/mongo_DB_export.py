from pymongo import MongoClient


def export(list_docs, collection_name, db_name):
    """
    Export data do a MongoDB server
    :param list_docs: list of dictionaries containing variant information
    :param collection_name: name of collection
    :param db_name: name of database
    :return: null
    """
    client = MongoClient()
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)
    collection.insert_many(list_docs)